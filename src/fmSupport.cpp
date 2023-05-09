#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "fmSupport.h"
#include "dy4.h"
#include <math.h>

using namespace std;

void fmDemod(std::vector<float> &I, std::vector<float> &Q, float &prevI, float &prevQ, std::vector<float> &fm_demod){
  //incomplete so far
	// prevI = 0.0
	// prevQ = 0.0

	// the default prev_phase phase is assumed to be zero, however
	// take note in block processing it must be explicitly controlled

		// empty vector to store the demodulated samples
		fm_demod.resize(I.size());

	  float currentI;
	  float currentQ;
	  float dI;
	  float dQ;

		// iterate through each of the I and Q pairs
		for (int k = 0; k < I.size(); k++){

			// the imaginary part (quadrature Q) and the real part (in-phase I)
			currentI = I[k];
			currentQ = Q[k];

			dI = currentI - prevI;
			dQ = currentQ - prevQ;

			// derivative of the phase
			if (pow(currentI,2) + pow(currentQ,2)== 0){
				fm_demod[k] = 0;
			} else {
				fm_demod[k] = (currentI*dQ - currentQ*dI)/(pow(currentI,2)+ pow(currentQ,2));
			}

			prevI = currentI;
			prevQ = currentQ;

	  }
}

void estimatePSD(std::vector<float> &samples, int Fs, std::vector<float> &freq, std::vector<float> &psd_est){

	int freq_bins = NFFT;
	float df = 1.0*Fs/freq_bins;

	freq.clear();
	psd_est.clear();

	for (float  i = 0; i < Fs/2; i+= df){

		freq.push_back(i);

	}


	std::vector<float> hann;

	for (int i = 0; i<freq_bins; i++){

		hann.push_back(std::pow(std::sin(i*PI/freq_bins),2));

	}


	std::vector<float> psd_list;

	int no_segments = (int)floor(samples.size()/(float)freq_bins);

	std::vector<float> windowed_samples;

	std::vector<complex<float>> Xf;

	std::vector<double> psd_seg;
	int j =0;

	for (int k = 0; k < no_segments; k++){

		j=0;
    windowed_samples.clear();
		for (int i = k*freq_bins; i < (k+1)*freq_bins; i++){

			windowed_samples.push_back(samples[i] * hann[j]);
			j++;

		}

		DFT(windowed_samples, Xf);

		Xf.resize((int)(freq_bins/2));

    psd_seg.clear();
		for (int i  = 0; i < (int)Xf.size(); i++){

			psd_seg.push_back((double)(2/((double)Fs*(double)freq_bins/2)) * (double) std::abs((std::pow(Xf[i],2))));

		}

		for (int i = 0; i < (int)psd_seg.size(); i++){

			psd_list.push_back((float)10*std::log10(psd_seg[i]));

		}

	}


	psd_est.resize((int)(freq_bins/2),0.0);

	for (int k  = 0; k < (int)(freq_bins/2); k++){

		for (int l = 0; l < no_segments; l++){

			psd_est[k] += psd_list[k+l*(int)(freq_bins/2)];

		}

		psd_est[k] = (psd_est[k]) / no_segments;

	}


}

void fmPLL(std::vector<float> &pllIn, std::vector<float> &ncoOut, std::vector<float> &state, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth){
	// scale factors for proportional/integrator terms
	// these scale factors were derived assuming the following:
	// damping factor of 0.707 (1 over square root of 2)
	// there is no oscillator gain and no phase detector gain
	float Cp = 2.666;
	float Ci = 3.555;

	// gain for the proportional term
	float Kp = normBandwidth*Cp;
	// gain for the integrator term
	float Ki = (normBandwidth*normBandwidth)*Ci;

	// output array for the NCO
	ncoOut.resize(pllIn.size(), 0.0);

	// initialize internal state
	float integrator = state[0];
	float phaseEst = state[1];
	float feedbackI = state[2];
	float feedbackQ = state[3];
	float trigOffset = state[5];
	ncoOut[0] = state[4];
	float trigArg;

	for (int k = 0; k < pllIn.size(); k++){

		// phase detector
		float errorI = pllIn[k]*1.0*feedbackI; // complex conjugate of the
		float errorQ = pllIn[k]*-1.0*feedbackQ; // feedback complex exponential

		// four-quadrnat arctangent discriminator for pahse error detection
		float errorD = std::atan2(errorQ,errorI);

		// loop filter
		integrator = integrator + Ki*errorD;
		// update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator;

		// internal oscillator
		trigOffset += 1.0;
		trigArg = 2.0*PI*(freq/Fs)*(trigOffset) + phaseEst;
		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);
		if (k+1 < ncoOut.size()){
				ncoOut[k+1] = std::cos(trigArg*ncoScale + phaseAdjust);
		}

		// for stereo only the in-phase NCO component should be returned
		// for block processing you should also return the state

	}


	// save state in vector
	state[0] = integrator;
	state[1] = phaseEst;
	state[2] = feedbackI;
	state[3] = feedbackQ;
	state[4] = std::cos(trigArg*ncoScale + phaseAdjust); // last value
	state[5] = trigOffset;

}

void plot_signals(std::vector<float> signal1, std::vector<float> signal2, std::vector<float> signal3, \
	std::vector<float> spectrum1, std::vector<float> spectrum2, std::vector<float> spectrum3){

		std::vector<float> freq;
		std::vector<float> timeVector;
		std::vector<float> psd_est;
		float T = 1.0/(240e3);

		for (float i = 0; i < signal3.size()*T; i+= T){

			timeVector.push_back(i);

		}

		//plot frequency domain signals
		estimatePSD(spectrum1, 240, freq, psd_est);
		logVector("spectrum1", freq, psd_est);

		estimatePSD(spectrum2, 240, freq, psd_est);
		logVector("spectrum2", freq, psd_est);

		estimatePSD(spectrum3, 240, freq, psd_est);
		logVector("spectrum3", freq, psd_est);

		//plot time domain signals
		logVector("signal1", timeVector, signal1);
		logVector("signal2", timeVector, signal2);
		logVector("signal3", timeVector, signal3);

		std::cerr << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/freq.gnuplot > ../data/freq.png\n";
		std::cerr << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/time.gnuplot > ../data/time.png\n";

}

//exactly the same as fmPLL but returns the quadrature component
void fmPLLQ(std::vector<float> &pllIn, std::vector<float> &ncoOut, std::vector<float> &state, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth){
	// scale factors for proportional/integrator terms
	// these scale factors were derived assuming the following:
	// damping factor of 0.707 (1 over square root of 2)
	// there is no oscillator gain and no phase detector gain
	float Cp = 2.666;
	float Ci = 3.555;

	// gain for the proportional term
	float Kp = normBandwidth*Cp;
	// gain for the integrator term
	float Ki = (normBandwidth*normBandwidth)*Ci;

	// output array for the NCO
	ncoOut.resize(pllIn.size(), 0.0);

	// initialize internal state
	float integrator = state[0];
	float phaseEst = state[1];
	float feedbackI = state[2];
	float feedbackQ = state[3];
	float trigOffset = state[5];
	ncoOut[0] = state[4];
	float trigArg;

	for (int k = 0; k < pllIn.size(); k++){

		// phase detector
		float errorI = pllIn[k]*1.0*feedbackI; // complex conjugate of the
		float errorQ = pllIn[k]*-1.0*feedbackQ; // feedback complex exponential

		// four-quadrnat arctangent discriminator for pahse error detection
		float errorD = std::atan2(errorQ,errorI);

		// loop filter
		integrator = integrator + Ki*errorD;
		// update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator;

		// internal oscillator
		trigOffset += 1.0;
		trigArg = 2.0*PI*(freq/Fs)*(trigOffset) + phaseEst;
		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);
		if (k+1 < ncoOut.size()){
				ncoOut[k+1] = std::sin(trigArg*ncoScale + phaseAdjust);
		}

		// for stereo only the in-phase NCO component should be returned
		// for block processing you should also return the state

	}


	// save state in vector
	state[0] = integrator;
	state[1] = phaseEst;
	state[2] = feedbackI;
	state[3] = feedbackQ;
	state[4] = std::sin(trigArg*ncoScale + phaseAdjust); // last value
	state[5] = trigOffset;

}
