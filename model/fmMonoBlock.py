#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD
# for take-home add your functions

rf_Fs = 2.4e6
rf_Fc = 100e3
audio_Fc = 16e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_Fs1 = 240e3
audio_decim = 5
# add other settings for audio, like filter taps, ...

# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs takehome (il_vs_th = 1)
il_vs_th = 1

def impulseResponse(Fs, Fc, N_taps):
	normCutoff = Fc/(Fs/2)
	h = np.zeros(N_taps)

	for i in range(N_taps -1):
		if i == (N_taps-1)/2:
			h[i] = normCutoff
		else:
			trash = math.pi*normCutoff*(i-((N_taps-1)/2))
			h[i] = normCutoff*math.sin(trash)/trash
		h[i] = h[i]*math.pow(math.sin((i*math.pi)/N_taps),2)
	return h

def convolution(h,x,zi=None):
	if len(x) < len(h):
		return

	if zi is None:
		y = np.zeros(len(x) + len(h) - 1)
		for n in range(len(y)):
			for k in range(len(h)):
				if (n-k>=0 and n-k<len(x)):
					y[n] = y[n] + x[n-k]*h[k]
		return y
	else:
		zf= x[-len(h):-1]
		y = np.zeros(len(x))
		for n in range(len(y)):
			for k in range(len(h)):
				if n-k >= 0:
					y[n] = y[n] + h[k]*x[n-k]
				else:
					y[n] = y[n] + h[k]*zi[n-k]
		return y, zf

#bruh
# def fmDemodTH(i_ds, q_ds, state_i,state_q, state_phase):
# 	h = np.array([1,0,-1])
# 	next_i = i_ds[-2:]
# 	next_q = q_ds[-2:]
# 	fm_demod = np.zeros(len(i_ds))
# 	denom = math.pow(i_ds[0],2)+math.pow(q_ds[0],2)
# 	if (denom != 0):
# 		fm_demod[0] = state_phase/denom
# 	for n in range(1,len(i_ds)):
# 		denom = math.pow(i_ds[n],2)+math.pow(q_ds[n],2)
# 		i_diff = 0
# 		q_diff = 0
# 		for k in range(3):
# 			if n-k >= 0:
# 				i_diff = i_diff + h[k]*i_ds[n-k]
# 				q_diff = q_diff + h[k]*q_ds[n-k]
# 			else:
# 				i_diff = i_diff + h[k]*state_i[n-k]
# 				q_diff = q_diff + h[k]*state_q[n-k]
#
#
# 		numerator = i_ds[n-1]*q_diff - q_ds[n-1]*i_diff
# 		if (denom != 0):
# 			fm_demod[n] = numerator/denom
#
# 	numerator = i_ds[-1]*q_diff - q_ds[-1]*i_diff

def fmDemodTH(i_ds, q_ds, state_i,state_q):
	next_i = i_ds[-1]
	next_q = q_ds[-1]
	fm_demod = np.zeros(len(i_ds))
	for n in range(0,len(i_ds)):
		denom = math.pow(i_ds[n],2)+math.pow(q_ds[n],2)
		if (n==0):
			i_diff = i_ds[n] - state_i
			q_diff = q_ds[n] - state_q
		else:
			i_diff = i_ds[n] - i_ds[n-1]
			q_diff = q_ds[n] - q_ds[n-1]

		numerator = i_ds[n]*q_diff - q_ds[n]*i_diff
		if (denom != 0):
			fm_demod[n] = numerator/denom

	return fm_demod, next_i, next_q

pll_state = {
        "integrator" : 0.0,
        "phaseEst" : 0.0,
        "feedbackI" : 1.0,
        "feedbackQ" : 0.0,
        "ncoOut_0" : 1.0,
        "trigOffset" : 0
}

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/iq_samples.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle  as for rf_coeff (but different arguments, of course)
		audio_coeff = signal.firwin(rf_taps, audio_Fc/(audio_Fs1/2), window=('hann'))
	else:
		# to be updated by you for the takehome exercise
		# with your own code for impulse response generation
		audio_coeff = impulseResponse(audio_Fs1, audio_Fc, rf_taps)

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_phase = 0
	# add state as needed for the mono channel filter
	state_i = 0
	state_q = 0

	# audio buffer that stores all the audio blocks
	audio_data = np.array([]) # used to concatenate filtered blocks (audio data)

	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < len(iq_data):

		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print('Processing block ' + str(block_count))

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
				zi=state_q_lpf_100k)

		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		if il_vs_th == 0:
			# already given to you for the in-lab
			# take particular notice of the "special" state-saving
			fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
		else:
			# you will need to implement your own FM demodulation based on:
			# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
			# see more comments on fmSupportLib.py - take particular notice that
			# you MUST have also "custom" state-saving for your own FM demodulator
			fm_demod, state_i, state_q = fmDemodTH(i_ds, q_ds, state_i, state_q)

		# extract the mono audio data through filtering
		if il_vs_th == 0:
			# to be updated by you during the in-lab session based on lfilter
			# same principle as for i_filt or q_filt (but different arguments)
			audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)
		else:
			# to be updated by you for the takehome exercise
			# with your own code for BLOCK convolution
			audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)

		# downsample audio data
		# to be updated by you during in-lab (same code for takehome)
		audio_block = audio_filt[::5]

		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data
		#
		audio_data = np.concatenate((audio_data, audio_block))
		#

		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
			ax0.clear()
			fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after extracting mono audio
			# ... change as needed

			# plot PSD of selected block after downsampling mono audio
			# ... change as needed

			# save figure to file
			fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = "../data/fmMonoBlock.wav"
	#for i in range(len(audio_data)):
		#print(audio_data[i])
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	# plt.show()
