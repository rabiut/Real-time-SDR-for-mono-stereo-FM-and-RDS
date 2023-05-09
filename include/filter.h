/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void impulseResponseBPF(float Fb, float Fe, float Fs, unsigned short int num_taps, std::vector<float> &h);
void impulseResponseAPF(float Fs, unsigned short int num_taps, std::vector<float> &h);
void APFnew(std::vector<float> &y, std::vector<float> &x, std::vector<float> &state);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void convolveFIR(std::vector<float> &y, std::vector<float> &x, std::vector<float> &h, std::vector<float> &zi);
void convolveFIR(std::vector<float> &y, std::vector<float> &x, std::vector<float> &h, std::vector<float> &zi, int d);
void convolveFIR(std::vector<float> &y, std::vector<float> &x, std::vector<float> &h, std::vector<float> &zi, int d, int u);

#endif // DY4_FILTER_H
