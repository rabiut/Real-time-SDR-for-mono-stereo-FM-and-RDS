/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <cmath>
#include <vector>


// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// bring your own functionality
  h.clear();
	h.resize(num_taps, 0.0); // allocate memory

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	float norm_cutoff = Fc/(Fs/2); // define normalized cutoff frequency

	for (int i = 0; i < num_taps; i++){
		if (i == (num_taps-1)/2){
			h[i] = norm_cutoff;
		} else {
			h[i] = norm_cutoff*(std::sin(PI*norm_cutoff*(i-(num_taps-1)/2))/(PI*norm_cutoff*(i-(num_taps-1)/2)));
		}
		h[i] = h[i]*std::pow(std::sin(i*PI/num_taps),2); // apply the Hann window n**2
	}
}

void impulseResponseBPF(float Fb, float Fe, float Fs, unsigned short int num_taps, std::vector<float> &h)
{
	// bring your own functionality
  h.clear();
	h.resize(num_taps, 0.0); // allocate memory

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	float norm_center = ((Fe+Fb)/2.0)/(Fs/2.0); // define normalized center frequency
  float norm_pass = (Fe-Fb)/(Fs/2.0); // define normalized pass frequency

	for (int i = 0; i < num_taps; i++){
		if (i == (num_taps-1.0)/2.0){
			h[i] = norm_pass;
		} else {
			h[i] = norm_pass*(std::sin(PI*norm_pass/2.0*(i-(num_taps-1.0)/2.0))/((PI*norm_pass/2.0)*(i-(num_taps-1)/2.0)));
		}
		h[i] = h[i]*std::cos(i*PI*norm_center); // apply a frequency shift by the center frequency
    h[i] = h[i]*std::pow(std::sin(i*PI/num_taps),2.0); // apply the Hann window n**2
	}
}

void impulseResponseAPF(float Fs, unsigned short int num_taps, std::vector<float> &h)
{

  //single impulse at num_taps
  h.clear();
	h.resize(num_taps, 0.0);
  h[h.size()-1] = 1.0;

}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// bring your own functionality
  y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	// h: impulse response (coeff), x: input

	for (int n = 0; n < y.size(); n++){
		for (int k = 0; k < h.size(); k++){
			if (n-k >= 0 && n-k < x.size()){
				y[n] += h[k] * x[n-k];
			}
		}
	}
}

//for block processing
void convolveFIR(std::vector<float> &y, std::vector<float> &x, std::vector<float> &h, std::vector<float> &zi)
{

		y.clear(); y.resize(x.size(), 0.0);
    for (int n = 0; n < y.size(); n++) {
        for (int k = 0; k < h.size(); k++){
            if (n - k >= 0) {
                y[n] += x[n - k] * h[k];
            } else {
                y[n] += zi[zi.size() + (n - k)] * h[k];
            }
        }
    }

    for (int i = x.size() - h.size(); i < x.size(); i++){

      zi[i - x.size() + h.size()] = x[i];

    }

}

//for block processing + downsampling
void convolveFIR(std::vector<float> &y, std::vector<float> &x, std::vector<float> &h, std::vector<float> &zi, int d)
{

		y.clear(); y.resize(x.size()/d, 0.0);
    int n1;
    for (int n = 0; n < x.size(); n += d) {
        n1 = n/d;
        for (int k = 0; k < h.size(); k++){
            if (n - k >= 0) {
                y[n1] += x[n - k] * h[k];
            } else {
                y[n1] += zi[zi.size() + (n - k)] * h[k];
            }
        }
    }

    for (int i = x.size() - h.size(); i < x.size(); i++){

      zi[i - x.size() + h.size()] = x[i];

    }

}

void convolveFIR(std::vector<float> &y, std::vector<float> &x, std::vector<float> &h, std::vector<float> &zi, int d, int u)
{

		y.clear(); y.resize((x.size()*u)/d, 0.0);
    int p;
    int i;
    for (int n = 0; n < y.size(); n++) {
        p = (n*d)%u;
        for (int k = p; k < h.size(); k+=u){
            i = (int)(n*d - k)/u;
            if (i >= 0) {
                y[n] += x[i] * h[k]*u;
            } else {
                y[n] += zi[zi.size() + i] * h[k]* u;
            }
        }
    }

    for (int i = x.size() - zi.size(); i < x.size(); i++){

      zi[i - x.size() + zi.size()] = x[i];

    }


}
