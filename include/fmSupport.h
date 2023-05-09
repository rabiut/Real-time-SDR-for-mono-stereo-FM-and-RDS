#ifndef DY4_FMSUPPORT_H
#define DY4_FMSUPPORT_H

// add headers as needed
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

void fmDemod(std::vector<float> &I, std::vector<float> &Q, float &prevI, float &prevQ, std::vector<float> &fm_demod);

void estimatePSD(std::vector<float> &samples, int Fs, std::vector<float> &freq, std::vector<float> &psd_est);

void fmPLL(std::vector<float> &pllIn, std::vector<float> &ncoOut, std::vector<float> &state, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth);

void plot_signals(std::vector<float> signal1, std::vector<float> signal2, std::vector<float> signal3, std::vector<float> spectrum1, std::vector<float> spectrum2, std::vector<float> spectrum3);

void fmPLLQ(std::vector<float> &pllIn, std::vector<float> &ncoOut, std::vector<float> &state, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth);

#endif // DY4_FOURIER_H
