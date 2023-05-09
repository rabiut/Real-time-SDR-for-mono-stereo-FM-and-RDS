/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "fmSupport.h"
#include "pthread.h"
#include "rdsSupport.h"
#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <algorithm>
#include <chrono>
#include <string>

// GLOBAL PARAMETERS
int mode;
int channel;
// initialize sampling rates and resampling factors, default: mode 0
float rf_Fs;
float rf_Fc;
float audio_Fc;
unsigned int rf_taps;
int rf_decim;

float IF_Fs;
int audio_decim;
int audio_up;
int mode2_fix;
int SPS;
int rds_up;
int rds_down;

int rds_symbol_rate;

unsigned int block_size;

//vectors that get saved and plotted
std::vector<float> spectrum1;
std::vector<float> spectrum2;
std::vector<float> spectrum3;
std::vector<float> signal1;
std::vector<float> signal2;
std::vector<float> signal3;
std::vector<float> scatter1;
std::vector<float> scatter2;
std::vector<float> bits1;


void rf_thread(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, std::condition_variable& my_cvar, std::atomic<int> &blocksPulled);
void audio_thread(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, std::condition_variable& my_cvar, std::atomic<int> &blocksPulled);
void rds_thread(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, std::condition_variable& my_cvar, std::atomic<int> &blocksPulled);


int main(int argc, char* argv[]) {

	mode = 0;
	channel = 0; // 0- mono, 1- stereo

	if (argc == 2){
		mode = atoi(argv[1]);
		if (mode > 3){
			std::cerr << "invalid mode: " << mode << std::endl;
			exit(1);
		}
	} else if (argc == 3){
		mode = atoi(argv[1]);
		if (mode > 3){
			std::cerr << "invalid mode: " << mode << std::endl;
			exit(1);
		}
		channel  = atoi(argv[2]);
		if (channel > 1){
			std::cerr << "invalid channel: " << channel << std::endl;
			exit(1);
		}
	} else if (argc < 2){
		std::cerr << "u wrong lol" << std::endl;
	}

	if (channel == 0){
		std::cerr << "Operating in mode "<< mode << " for mono" <<std::endl;
	} else if (channel == 1){
		std::cerr << "Operating in mode "<< mode << " for stereo" <<std::endl;
	}

	// MODE PARAMETERS
	rf_Fs = 2.4e6;
	rf_Fc = 100e3;
	audio_Fc = 16e3;
	rf_taps = 151;
	rf_decim = 10;

	IF_Fs = 240e3;
	audio_decim = 5;
	audio_up = 1;
	mode2_fix = 1;
	SPS = 28;
	rds_up = 133;
	rds_down = 480;

	// specify for modes 1-3
	if (mode == 1) {
		rf_Fs = 960e3;
		rf_decim = 4;
		audio_decim = 5;

	} else if (mode == 2) {
		audio_up = 147;
		audio_decim = 800;
		mode2_fix = 6;
		SPS = 45;
		rds_up = 57;
		rds_down = 128;

	} else if (mode == 3) {
		rf_Fs = 1.152e6;
		IF_Fs = 384e3;
		rf_decim = 3;
		audio_up = 147;
		audio_decim = 1280;

	}

	rds_symbol_rate = 2375*SPS;


	// determine block processing size
	// select a block_size that is a multiple of KB
	// and a multiple of decimation factors
	block_size = (1024 * rf_decim * audio_decim * 2* mode2_fix)/audio_up;

	std::queue<std::vector<float>> my_queue;
	std::mutex my_mutex;
	std::condition_variable my_cvar;
	std::atomic<int> blocksPulled;

	//------------------------------ BEGIN PROCESSING (THREADING) - READING DATA STREAM--------------------------------------------//
	// rf_thread
	std::thread RF_producer_thread = std::thread(rf_thread, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), std::ref(blocksPulled));
	// audio_thread
	std::thread AUDIO_consumer_thread = std::thread(audio_thread, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), std::ref(blocksPulled));
	// rds thread
	std::thread RDS_consumer_thread = std::thread(rds_thread, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), std::ref(blocksPulled));

	RF_producer_thread.join();
	AUDIO_consumer_thread.join();
	RDS_consumer_thread.join();

	//---------------------MISC ----------------------//

	return 0;

} /***********END OF MAIN************/


/**************************************************************PRODUCER: rf_thread()********************************************************************/
void rf_thread(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, std::condition_variable& my_cvar, std::atomic<int>& blocksPulled){

	std::cerr << "Beginning RF Thread..." << std::endl;
	//---------------------RF DEMODULATION COMPONENTS--------------------------//
	// coefficients for the low-pass rf filter
	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);

	// i & q states needed for continuity in block processing
	std::vector<float>state_i_rf;
	std::vector<float>state_q_rf;
	state_i_rf.clear(); state_i_rf.resize(rf_taps-1, 0.0);
	state_q_rf.clear(); state_q_rf.resize(rf_taps-1, 0.0);

	//iq samples to demodulated vectors
	float state_i = 0.0;
	float state_q = 0.0;
	std::vector<float> i_block;
	std::vector<float> q_block;
	std::vector<float> fm_demod;
	std::vector<float> i_filt;
	std::vector<float> q_filt;

	blocksPulled = 0;

	std::cerr << "Starting Audio Processing" << std::endl;

	for (unsigned int block_id = 0; ;block_id++) {
		//-------------- BEGIN RF DEMODULATION-------------//
		std::vector<float> block_data(block_size);
		readStdinBlockData(block_size, block_id, block_data);

		//stop when there is no more data
		if ((std::cin.rdstate()) != 0){
			std::cerr << "End of input steam" << std::endl;
			std::vector<float> sampleVector;
			for (int i = 0; i < bits1.size(); i++){

				sampleVector.push_back(i);

			}
			//phase tuning plots
			logVector("IQ", scatter1, scatter2);
			logVector("symbPairing", sampleVector,bits1);
			plot_signals(signal1,signal2,signal3,spectrum1,spectrum2,spectrum3);
			exit(1);
		}

		//split data into i and q samples
		i_block.clear();
		q_block.clear();
		for (int i = 0; i < block_data.size(); i++){
			if (i%2 == 0){
				i_block.push_back(block_data[i]);
			} else{
				q_block.push_back(block_data[i]);
			}
		}

		//first LPF to 100k
		convolveFIR(i_filt, i_block, rf_coeff, state_i_rf, rf_decim);
		convolveFIR(q_filt, q_block, rf_coeff, state_q_rf, rf_decim);

		//fm_demodulation
		fmDemod(i_filt, q_filt, state_i, state_q, fm_demod);

		//-------------- END RF DEMODULATION-------------//

		std::unique_lock<std::mutex> my_lock(my_mutex);
		while (my_queue.size()>=QUEUE_ELEMS || blocksPulled != 0) { // if the queue is full
			my_cvar.wait(my_lock); // wait
		}

		//change elem to vector of iq samples
		my_queue.push(fm_demod);
		blocksPulled++;
		//std::cerr << "Producer pushes to queue..." << std::endl;
		my_cvar.notify_all();
		my_lock.unlock(); // unlock

	}
}

/**************************************************************CONSUMER: audio_thread()********************************************************************/
void audio_thread(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, std::condition_variable& my_cvar, std::atomic<int>& blocksPulled){


	std::cerr << "Beginning Audio Thread..." << std::endl;
	//---------------------AUDIO COMPONENTS----------------------------//
	//lpf to extract mono
	std::vector<float> audio_coeff;
	impulseResponseLPF(IF_Fs*audio_up, audio_Fc, rf_taps*audio_up*2, audio_coeff);

	// coefficients for the band-pass filter (pilot & audio channel)
	std::vector<float> stereo_pilot_coeff;
	impulseResponseBPF(18.5e3, 19.5e3, IF_Fs, rf_taps, stereo_pilot_coeff);
	std::vector<float> stereo_channel_coeff;
	impulseResponseBPF(22e3, 54e3, IF_Fs, rf_taps, stereo_channel_coeff);

	//LPF after demodulating stereo channel
	std::vector<float> DSBSC_coeff;
	impulseResponseLPF(IF_Fs*audio_up, 16e3, rf_taps*audio_up, DSBSC_coeff);
	std::vector<float> DSBSC_filt_state;
	DSBSC_filt_state.clear(); DSBSC_filt_state.resize(rf_taps-1, 0.0);

	// get mono filter state
	std::vector<float>monoFilterState;
  monoFilterState.clear(); monoFilterState.resize((rf_taps*2)-1, 0.0);

	// get stereo pilot filter state
	std::vector<float>SPilotFilterState;
  SPilotFilterState.clear(); SPilotFilterState.resize(rf_taps-1, 0.0);

	// get stereo channel filter state
	std::vector<float>SChannelFilterState;
  SChannelFilterState.clear(); SChannelFilterState.resize(rf_taps-1, 0.0);

	//PLL states plus output
	std::vector<float> clean_stereo_pilot;
	std::vector<float> pll_state;
	pll_state.clear(); pll_state.resize(6,0.0);
	pll_state[2] = 1.0;
	pll_state[4] = 1.0;

	//intermediate processed stereo data
	std::vector<float> stereo_signal, stereo_pilot, stereo_channel;

	//final audio data
	std::vector<float> processed_data;
	std::vector<float> stereo_filtered;
	std::vector<float> Left, Right;
	std::vector<short int> audio_data;

	//CONTINUE TO EXTRACT EACH BLOCK FROM QUEUE SO LONG THERE IS SOMETHING TO "CONSUME"
	while (1) {
		std::unique_lock<std::mutex> my_lock(my_mutex);
		while (my_queue.empty() || blocksPulled != 1) {
			my_cvar.wait(my_lock);
		}

		//extract block from queue
		std::vector<float> fm_demod = my_queue.front();
		//my_queue.pop();
		blocksPulled++;
		my_cvar.notify_all();
		my_lock.unlock();
		//std::cerr << "Consumer audio extracts from queue" << std::endl;

		//---------------------BEGIN AUDIO (MONO & STEREO) PROCESSING ---------------------//
		//MONO audio extraction
		convolveFIR(processed_data, fm_demod, audio_coeff, monoFilterState, audio_decim, audio_up);
		//STEREO pilot extraction
		convolveFIR(stereo_pilot, fm_demod, stereo_pilot_coeff, SPilotFilterState);
		//stereo channel extraction
		convolveFIR(stereo_channel, fm_demod, stereo_channel_coeff, SChannelFilterState);
		//pll
		fmPLL(stereo_pilot, clean_stereo_pilot, pll_state, 19e3, IF_Fs, 2.0, 0.0, 0.01);

		stereo_signal.clear();

		for (int i = 0; i <  clean_stereo_pilot.size(); i++){
			stereo_signal.push_back(2*clean_stereo_pilot[i]*stereo_channel[i]);
		}

		convolveFIR(stereo_filtered, stereo_signal, DSBSC_coeff, DSBSC_filt_state, audio_decim, audio_up);

		//L channel processed_data + stereo_filtered
		Left.clear(); Right.clear();
		for (int i = 0; i <  processed_data.size(); i++){
				Left.push_back(processed_data[i] + stereo_filtered[i]);
				Right.push_back(stereo_filtered[i] - processed_data[i]);
		}

		//---------------------END AUDIO (MONO & STEREO) PROCESSING -------------------------//

		//---------------------BEGIN AUDIO WRITING -------------------------//
		audio_data.clear();

		if (channel == 0){
			// mono writing
			for (unsigned int k =0; k < processed_data.size(); k++){
				if (std::isnan(processed_data[k])){
					audio_data.push_back(0);
				} else {
					audio_data.push_back(static_cast<short int>(processed_data[k]*16384));
				}
			}
		} else if (channel == 1){
			//stereo writing
			for (unsigned int k =0; k < Left.size()*2; k++){
				if (k%2 == 0){
					if (std::isnan(Right[k/2])){
						audio_data.push_back(0);
					} else {
						audio_data.push_back(static_cast<short int>(Right[k/2]*16384));
					}
				} else {
					if (std::isnan(Left[(k-1)/2])){
						audio_data.push_back(0);
					} else {
						audio_data.push_back(static_cast<short int>(Left[(k-1)/2]*16384));
					}
				}
			}
		}

		//write the block of data
		fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);

	}

}

/**************************************************************CONSUMER: rds_thread()********************************************************************/
void rds_thread(std::queue<std::vector<float>> &my_queue, std::mutex& my_mutex, std::condition_variable& my_cvar, std::atomic<int>& blocksPulled){

	// -----------------RDS STUFF -----------------------//

	int block_id = 0;

	//BPF for RDS
	std::vector<float> rds_coeff;
	impulseResponseBPF(54e3, 60e3, IF_Fs, rf_taps, rds_coeff);
	std::vector<float> rds_bpf_state;
	rds_bpf_state.clear(); rds_bpf_state.resize(rf_taps-1, 0.0);

	//APF for delay
	std::vector<float> apf_coeff;
	impulseResponseAPF(IF_Fs, (rf_taps-1)/2, apf_coeff);
	std::vector<float> apf_state;
	apf_state.clear(); apf_state.resize((rf_taps-1)/2 - 1, 0.0);

	//BPF for RDS carrier
	std::vector<float> rds_carrier_coeff;
	impulseResponseBPF(113.5e3, 114.5e3, IF_Fs, rf_taps, rds_carrier_coeff);
	std::vector<float> rds_carrier_bpf_state;
	rds_carrier_bpf_state.clear(); rds_carrier_bpf_state.resize(rf_taps-1, 0.0);

	//LPF after demod
	std::vector<float> rds_demod_coeff;
	impulseResponseLPF(IF_Fs*rds_up, 3e3, rf_taps*rds_up, rds_demod_coeff);
	std::vector<float> rds_demod_state;
	rds_demod_state.clear(); rds_demod_state.resize(rf_taps-1, 0.0);

	//pll stuff for carrier
	std::vector<float> rds_pll_state;
	rds_pll_state.clear(); rds_pll_state.resize(6,0.0);
	rds_pll_state[2] = 1.0;
	rds_pll_state[4] = 1.0;

	//roor raised cosine filter
	std::vector<float> rrc_coeffs;
	impulseResponseRootRaisedCosine(rds_symbol_rate, rf_taps, rrc_coeffs);
	std::vector<float> rrc_state;
	rrc_state.clear(); rrc_state.resize(rf_taps-1,0.0);

	//stores different stages of processed rds data
	std::vector<float> rds_raw;
	std::vector<float> rds_squared;
	std::vector<float> rds_carrier;
	std::vector<float> rds_carrier_clean;
	std::vector<float> rds_delayed;
	std::vector<float> rds_demod;
	std::vector<float> rds_message;
	std::vector<float> rrc_filtered;
	std::vector<bool> bitstreamRDS;
	std::vector<bool> largeBitsream;
	std::vector<bool> windowState;
	std::vector<float> fm_demod_large;
	std::string PService = "";
	int DIprev = -1;
	int DI = 0;
	bool bitState = 0;
	int symbolState = -1;
	bool synced = false;
	int symbCount = 0;

	//used just for phase matching(not actually necessary for final product)
	std::vector<float> rds_carrier_cleanQ;
	std::vector<float> rds_demodQ;
	std::vector<float> rds_messageQ;
	std::vector<float> rrc_filteredQ;
	std::vector<float> rds_pll_stateQ;
	rds_pll_stateQ.clear(); rds_pll_stateQ.resize(6,0.0);
	rds_pll_stateQ[2] = 1.0;
	rds_pll_stateQ[4] = 1.0;
	std::vector<float> rrc_stateQ;
	rrc_stateQ.clear(); rrc_stateQ.resize(rf_taps-1,0.0);
	std::vector<float> rds_demod_stateQ;
	rds_demod_stateQ.clear(); rds_demod_stateQ.resize(rf_taps-1, 0.0);


	//this waits for unlock i think
	while (1) {
		std::unique_lock<std::mutex> my_lock(my_mutex);
		while (my_queue.empty() || blocksPulled != 2) {
			my_cvar.wait(my_lock);
		}

		//replace to take in iq samples
		std::vector<float> fm_demod = my_queue.front();
		blocksPulled=0;
		my_queue.pop();
		my_cvar.notify_all();
		my_lock.unlock();

		//std::cerr << "Consumer rds extracts from queue"<< std::endl;

		//---------------------BEGIN RDS PROCESSING ---------------------//
		for (int i = 0; i < fm_demod.size(); i++){fm_demod_large.push_back(fm_demod[i]);}
		if (block_id % 5 == 0){

			//rds cahnnel extraction
			convolveFIR(rds_raw, fm_demod_large, rds_coeff, rds_bpf_state);
			fm_demod_large.clear();

			//RDS delay just so it matches the delay introduced during carrier extraction
			convolveFIR(rds_delayed, rds_raw, apf_coeff, apf_state);
			//APFnew(rds_delayed, rds_raw, apf_state);

			//squaring nonlinearity
			rds_squared.clear();
			for (int i = 0; i < rds_raw.size(); i++){

				rds_squared.push_back(rds_raw[i]*rds_raw[i]);

			}

			//carrier extraction
			convolveFIR(rds_carrier, rds_squared, rds_carrier_coeff, rds_carrier_bpf_state);

			//cleaning the carrier signal
			fmPLL(rds_carrier, rds_carrier_clean, rds_pll_state, 114e3, IF_Fs, 0.5, -(1.1*PI)/8.0, 0.001);

			//quadrature
			fmPLLQ(rds_carrier, rds_carrier_cleanQ, rds_pll_stateQ, 114e3, IF_Fs, 0.5, -(1.1*PI)/8.0, 0.001);

			//mixer
			rds_demod.clear();
			rds_demodQ.clear();
			for (int i = 0; i <rds_carrier_clean.size(); i++){

				rds_demod.push_back(2.0*rds_carrier_clean[i]*rds_delayed[i]);
				rds_demodQ.push_back(2.0*rds_carrier_cleanQ[i]*rds_delayed[i]);

			}

			//LPF to get rid of rds image
			convolveFIR(rds_message, rds_demod, rds_demod_coeff, rds_demod_state, rds_down, rds_up);
			//quadrature
			convolveFIR(rds_messageQ, rds_demodQ, rds_demod_coeff, rds_demod_stateQ, rds_down, rds_up);

			//root raised cosine filter
			convolveFIR(rrc_filtered, rds_message, rrc_coeffs, rrc_state);
			//quadrature
			convolveFIR(rrc_filteredQ, rds_messageQ, rrc_coeffs, rrc_stateQ);

			//finds start sample every time cuz a predefined one would get worse over time
			findPeakSample(rrc_filtered, symbCount, SPS);

			//lots of samples and blocks are needed for plotting constellation
			IQSampler(rrc_filtered, rrc_filteredQ, scatter1, scatter2, symbCount, SPS);

			//if its taking too long to sync then it started on the wrong symbol
			//so it just removes the first symbol to re-align
			if (!synced && block_id % 15 == 0){
				for(int i=0;i<SPS;i++){rrc_filtered.erase(rrc_filtered.begin());}
			}

			//Gets symbols, manchestor decodes, and differntial decodes
			differntialDecode(rrc_filtered, bitstreamRDS, symbolState, bitState, symbCount, SPS);

			//sees if it can syncronize to any frame within bitStreamRDS
			frameSync(bitstreamRDS, windowState, synced, PService, DIprev, DI);

			//saving data for plotting
			if (block_id == 50){
				spectrum1 = rds_raw;
				spectrum2 = rds_carrier;
				spectrum3 = rds_demod;

				signal1 = rrc_filteredQ;
				signal2 = rds_message;
				signal3 = rrc_filtered;
				//exit(1);

			}

			//std::cerr << "processed rds for block: " << block_id/5 << std::endl;

		}
		block_id++;

		//---------------------END RDS PROCESSING ---------------------//

	}

}
