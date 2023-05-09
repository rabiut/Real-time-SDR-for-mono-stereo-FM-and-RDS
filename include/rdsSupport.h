#ifndef DY4_RDSSUPPORT_H
#define DY4_RDSSUPPORT_H

// add headers as needed
#include <iostream>
#include <vector>
#include <string>

void impulseResponseRootRaisedCosine(float Fs, int n_taps, std::vector<float> &impulseResponseRRC);

void differntialDecode(std::vector<float> &preCDR, std::vector<bool> &decodedBitstream, int &symbolState, bool &bitState, int symbCount, int SPS);

void checkFrame(std::vector<bool> &bitStream, int &blockType);

void frameSync(std::vector<bool> &bitStream, std::vector<bool> &windowState, bool &synced, std::string &PService, int &DIprev, int &DI);

void findPeakSample(std::vector<float> &rrc_filtered, int &symbCount, int SPS);

void IQSampler(std::vector<float> &preCDRI, std::vector<float> &preCDRQ,std::vector<float> &Isamples, std::vector<float> &Qsamples, int symbCount, int SPS);

void binarytoString(std::vector<bool> &bits, int blockType);

char binToChar(std::vector<bool> binary);

#endif // DY4_FOURIER_H
