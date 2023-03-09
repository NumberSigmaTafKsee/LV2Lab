#pragma once 
#include "Undenormal.hpp"
#include <cmath>

namespace Analog::Filters::StateVariableFilter2
{
    struct StateVariableFilter : public FilterProcessor
    {
        /*
        cutoff = cutoff freq in Hz
        fs = sampling frequency //(e.g. 44100Hz)
        f = 2 sin (pi * cutoff / fs) //[approximately]
        q = resonance/bandwidth [0 < q <= 1]  most res: q=1, less: q=0
        low = lowpass output
        high = highpass output
        band = bandpass output
        notch = notch output
        */
        DspFloatType cutoff,scale,fs,low,high,band,notch;
            
        StateVariableFilter(DspFloatType Fc, DspFloatType Q, DspFloatType Fs) : FilterProcessor() {
            scale = Q;
            cutoff= Fc;
            fs    = Fs;
            low=high=band=notch=0;
        }
        void setCutoff(DspFloatType F) { cutoff = F; }
        void setResonance(DspFloatType R) { scale = 1.25*(1.0-R); }

        enum {
            PORT_CUTOFF,
            PORT_RESONANCE
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_RESONANCE: setResonance(v); break;
            }
        }        
        DspFloatType Tick(DspFloatType I, DspFloatType A = 1, DspFloatType X=0, DspFloatType Y=0)
        {
            Undenormal denormal;
            DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
            //--beginloop
            //I = tanhify(I);
            low = low + f * band;
            high = scale * I - low - scale*band;
            band = f * high + band;
            notch = high + low;
            return low;
        }
        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out) {
            Undenormal denormal;
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++) {
                DspFloatType f     = 2 * std::sin(2 * M_PI * cutoff/fs);        
                //--beginloop
                //I = tanhify(I);
                low = low + f * band;
                high = scale * in[i] - low - scale*band;
                band = f * high + band;
                notch = high + low;
                out[i] = low;
             }
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            ProcessSIMD(n,in,out);
        }
        void ProcessInplace(size_t n, DspFloatType * out) {
            ProcessSIMD(n,out,out);
        }
    };
}
