#include "dy4.h"
#include "rdsSupport.h"
#include <math.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <string>

bool p[26][10] = {{1,0,0,0,0,0,0,0,0,0}, \
    {0,1,0,0,0,0,0,0,0,0}, \
    {0,0,1,0,0,0,0,0,0,0}, \
    {0,0,0,1,0,0,0,0,0,0}, \
    {0,0,0,0,1,0,0,0,0,0}, \
    {0,0,0,0,0,1,0,0,0,0}, \
    {0,0,0,0,0,0,1,0,0,0}, \
    {0,0,0,0,0,0,0,1,0,0}, \
    {0,0,0,0,0,0,0,0,1,0}, \
    {0,0,0,0,0,0,0,0,0,1}, \
    {1,0,1,1,0,1,1,1,0,0}, \
    {0,1,0,1,1,0,1,1,1,0}, \
    {0,0,1,0,1,1,0,1,1,1}, \
    {1,0,1,0,0,0,0,1,1,1}, \
    {1,1,1,0,0,1,1,1,1,1}, \
    {1,1,0,0,0,1,0,0,1,1}, \
    {1,1,0,1,0,1,0,1,0,1}, \
    {1,1,0,1,1,1,0,1,1,0}, \
    {0,1,1,0,1,1,1,0,1,1}, \
    {1,0,0,0,0,0,0,0,0,1}, \
    {1,1,1,1,0,1,1,1,0,0}, \
    {0,1,1,1,1,0,1,1,1,0}, \
    {0,0,1,1,1,1,0,1,1,1}, \
    {1,0,1,0,1,0,0,1,1,1}, \
    {1,1,1,0,0,0,1,1,1,1}, \
    {1,1,0,0,0,1,1,0,1,1}};

std::vector<bool> A{1,1,1,1,0,1,1,0,0,0};
std::vector<bool> B{1,1,1,1,0,1,0,1,0,0};
std::vector<bool> C{1,0,0,1,0,1,1,1,0,0};
std::vector<bool> Cp{1,1,1,1,0,0,1,1,0,0};
std::vector<bool> D{1,0,0,1,0,1,1,0,0,0};

std::vector<bool> test{1,1,0,0, 0,0,1,0, 0,1,1,0, 0,0,0,0, 1,1,1,1, 0,1,1,0, 0,0};

void impulseResponseRootRaisedCosine(float Fs, int n_taps, std::vector<float> &impulseResponseRRC){

  float T_symbol = 1.0/2375.0;
  float beta = 0.9;
  impulseResponseRRC.clear();
  impulseResponseRRC.resize(n_taps,0.0);

  for (int k = 0; k < n_taps; k++){

    float t = (k-(n_taps/2.0))/Fs;
    if (t == 0.0){
      impulseResponseRRC[k] = 1.0 + beta*((4.0/PI)-1.0);
    } else if((t == -T_symbol/(4.0*beta)) || (t == T_symbol/(4.0*beta))) {
      impulseResponseRRC[k] = (beta/std::sqrt(2.0))*( (1.0+2.0/PI) * (std::sin(PI/(4.0*beta))) + ((1.0-2.0/PI) * std::cos(PI/(4.0*beta)) ));
    } else {
      impulseResponseRRC[k] = (std::sin(PI*t*(1.0-beta)/T_symbol) +  4.0*beta*(t/T_symbol)*std::cos(PI*t*(1.0+beta)/T_symbol)) /
      (PI*t*(1.0-(4.0*beta*t/T_symbol)*(4.0*beta*t/T_symbol))/T_symbol);
    }

  }


}

void differntialDecode(std::vector<float> &preCDR, std::vector<bool> &decodedBitstream, int &symbolState, bool &bitState, int symbCount, int SPS){

  decodedBitstream.clear();
  std::vector<bool> symbols;
  std::vector<bool> manchesterDecode;

  if (symbolState != -1){
    symbols.push_back((bool)symbolState);
  }

  //extracting symbols from post RRC filtered signal
  //symbcount is used to sample every SPS samples and maintian proper sampling between blocks
  while (symbCount < preCDR.size()){

    if (preCDR[symbCount] > 0){
      symbols.push_back(1);
    } else {
      symbols.push_back(0);
    }
    symbCount += SPS;

  }

  if (symbols.size()%2 != 0){
    symbolState = (bool)symbols[symbols.size()-1];
  } else {
    symbolState = -1;
  }

  //manchester decoding the symbols
  for (int i = 0; i < symbols.size() - 1; i+= 2){

      //pairing based off symbol transitions
      if (symbols[i] == 1){
        manchesterDecode.push_back(1);
      } else {
        manchesterDecode.push_back(0);
      }

  }

  //differntial decoding by XOR
  //last bit of manchester decoded bitstream is state saved for next block
  decodedBitstream.push_back(bitState ^ manchesterDecode[0]);
  for (int i = 1; i < manchesterDecode.size(); i++){
    decodedBitstream.push_back(manchesterDecode[i] ^ manchesterDecode[i-1]);
  }
  bitState = manchesterDecode[manchesterDecode.size() - 1];

}

void checkFrame(std::vector<bool> &bitStream, int &blockType){

  std::vector<bool> checkWord; checkWord.resize(10,0.0);

  //matrix multiplication
  for (int i = 0; i < 10; i++){

    for (int j  = 0; j < 26; j++){

      checkWord[i]  = checkWord[i] ^ (bitStream[j]*p[j][i]);

    }

  }

  //checking which block
  blockType = 0;
  if (checkWord == A){
    blockType = 1;
  } else if (checkWord == B){
    blockType = 2;
  } else if (checkWord == C){
    blockType = 3;
  } else if (checkWord == Cp){
    blockType = 4;
  } else if (checkWord == D){
    blockType = 5;
  }

}

void frameSync(std::vector<bool> &bitStream, std::vector<bool> &windowState, bool &synced, std::string &PService, int &DIprev, int &DI){

  int blockPos = 0;
  int blockType = 0;
  std::vector<bool> largeBitsream;
  std::vector<bool> window;
  std::vector<bool> c1;
  std::vector<bool> c2;

  //combines previous block values with input bitsream
  //the window already starts at the proper location that needs to be checked
  for (int i = 0; i < windowState.size(); i++){
    largeBitsream.push_back(windowState[i]);
  }
  for (int i = 0; i < bitStream.size(); i++){
    largeBitsream.push_back(bitStream[i]);
  }

  //keep checking windows as long as the bitstream allows
  while (blockPos + 26 < largeBitsream.size()){

      //prepares the window to be checked
      window.clear();
      for (int i = blockPos; i < blockPos + 26; i++){
        window.push_back(largeBitsream[i]);
      }
      //checks the block type based of the parity matrix
      checkFrame(window, blockType);

      //sync failed
      if (blockType == 0){
        blockPos++;
        //synced = false;
      } else {

        //sync succeeded and shifts window by 26 now that it is syncronized
        //std::cerr << "SYNCED! block: " << blockType << "  ";
        synced = true;
        blockPos += 26;
        if (blockType == 1){
          std::cerr << "PI Code: ";
          binarytoString(window, blockType);
          std::cerr << std::endl;

        } else if (blockType == 2){

          std::cerr << "Program Type: ";
          binarytoString(window, blockType);
          std::cerr << std::endl;

          //binary to decimal
          DI = int(window[14])*2 + int(window[15]);
/*
          std::cerr << " Group Type: ";
          for (int i = 0; i < 5; i ++){
            std::cerr << window[i];
          }
*/
        } else if (blockType == 5){

          //if we have the next 2 characters then read the block D data
          if (DI-1 == DIprev){
            c1.clear();
            c2.clear();

            for (int i = 0; i < 8; i ++){
              c1.push_back(window[i]);
            }

            for (int i = 8; i < 16; i ++){
              c2.push_back(window[i]);
            }

            PService += binToChar(c1);
            PService += binToChar(c2);
            //one message is constructed, output the string
            if (DI == 3){
              std::cerr << "Program Service: " << PService << std::endl;
            }
          }

          DIprev = DI;

        }

      }

    }

  //state saving last unused elements
  windowState.clear();
  for (int i = blockPos; i < largeBitsream.size(); i++){
      windowState.push_back(largeBitsream[i]);
  }


}

void findPeakSample(std::vector<float> &rrc_filtered, int &symbCount, int SPS){

  float totalDist = 0;
  std::vector<float> distances;
  std::vector<float> symbol;
  int numElems = 0;
  //goes through all possible start points
  for (int i = 0; i < SPS; i ++){

    totalDist = 0.0;
    numElems = 0;
    //gets the total distance of each sample from the origin
    for (int j = i; j < rrc_filtered.size(); j+=SPS){
      totalDist += abs(rrc_filtered[j]);
      numElems++;
    }

    distances.push_back(totalDist/numElems);

  }

  //gets the max distance because it would mean that the samples have the largest distances from 0
  float max = distances[0];
  for (int i  = 1; i < distances.size(); i++){
    if (distances[i] > max){
      max = distances[i];
      symbCount = i;
    }
  }


}

//only used for phase tuning/constellation plot data extraction
void IQSampler(std::vector<float> &preCDRI, std::vector<float> &preCDRQ,std::vector<float> &Isamples, std::vector<float> &Qsamples, int symbCount, int SPS){

  //extracting symbols from post RRC filtered signal
  //symbcount is used to sample every SPS samples and maintian proper sampling between blocks
  while (symbCount < preCDRI.size()){
    Isamples.push_back(preCDRI[symbCount]);
    Qsamples.push_back(preCDRQ[symbCount]);
    symbCount += SPS;
  }

}

void binarytoString(std::vector<bool> &bits, int blockType){
  std::string pI = ""; // PI code
  std::string PTY = ""; // Program Type

  for (int i = 0; i < 16; i ++){
    if (bits[i]) {
        pI += "1";
    } else {
        pI += "0";
    }
  }

  for (int i = 6; i < 11; i ++){
    if (bits[i]) {
        PTY += "1";
    } else {
        PTY += "0";
    }
  }

  // Program Type
  if (PTY =="00000"){ //0
    PTY = "Undefined";
  } else if (PTY == "00001"){ //1
    PTY = "News";
  } else if (PTY == "00010"){ //2
    PTY = "Information";
  } else if (PTY == "00011"){ //3
    PTY = "Sports";
  } else if (PTY == "00100"){ //4
    PTY = "Talk";
  } else if (PTY == "00101"){ //5
    PTY = "Rock";
  } else if (PTY == "00110"){ //6
    PTY = "Classic Rock";
  } else if (PTY == "00111"){ //7
    PTY = "Adult Hits";
  } else if (PTY == "01000"){ //8
    PTY = "Soft Rock";
  } else if (PTY == "01001"){ //9
    PTY = "Top 40";
  } else if (PTY == "01010"){ //10
    PTY = "Country";
  } else if (PTY == "01011"){ //11
    PTY = "Oldies";
  } else if (PTY == "01100"){ //12
    PTY = "Soft";
  } else if (PTY == "01101"){ //13
    PTY = "Nostalgia";
  } else if (PTY == "01110"){ //14
    PTY = "Jazz";
  } else if (PTY == "01111"){ //15
    PTY = "Classical";
  } else if (PTY == "10000"){ //16
    PTY = "Rhythm & Blues";
  } else if (PTY == "10001"){ //17
    PTY = "Soft Rhythm & Blues";
  } else if (PTY == "10010"){ //18
    PTY = "Foreign Language";
  } else if (PTY == "10011"){ //19
    PTY = "Religious Music";
  } else if (PTY == "10100"){ //20
    PTY = "Religious Talk";
  } else if (PTY == "10101"){ //21
    PTY = "Personality";
  } else if (PTY == "10110"){ //22
    PTY = "Public";
  } else if (PTY == "10111"){ //23
    PTY = "College";
  } else if (PTY == "11101"){ //29
    PTY = "Weather";
  } else if (PTY == "11110"){ //30
    PTY = "Emergency Test";
  } else if (PTY == "11111"){ //31
    PTY = "Emergency";
  } else { //24-28
    PTY = "Unassigned";
  }

  if (blockType == 1){
    pI = std::string(pI.length() % 4 ? 4 - pI.length() % 4 : 0, '0') + pI;
    std::unordered_map<std::string, char> hex_dict = {
        {"0000", '0'}, {"0001", '1'}, {"0010", '2'}, {"0011", '3'},
        {"0100", '4'}, {"0101", '5'}, {"0110", '6'}, {"0111", '7'},
        {"1000", '8'}, {"1001", '9'}, {"1010", 'a'}, {"1011", 'b'},
        {"1100", 'c'}, {"1101", 'd'}, {"1110", 'e'}, {"1111", 'f'}
    };
    std::string pIhex;
    for (size_t i = 0; i < pI.length(); i += 4) {
        std::string group = pI.substr(i, 4);
        pIhex += hex_dict[group];
    }
    std ::cerr << pIhex;

  } else if (blockType == 2){
    std::cerr << PTY;
  }


}

char binToChar(std::vector<bool> binary){

  int decimal = 0;
  for (int i = 0; i < binary.size(); i++){
    decimal += std::pow(2,i)*int(binary[7-i]);
  }

  return char(decimal);

}
