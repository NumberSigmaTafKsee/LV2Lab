#pragma once

#include "SoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"
#include "VABlitSawOscillator.hpp"

namespace Analog::Oscillators
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // square is made from subtracting out of phase sawtooth waves
    ///////////////////////////////////////////////////////////////////////////////////////
    struct BlitSquare : public OscillatorProcessor
    {
        FX::Filters::OnePole block;
        Analog::Oscillators::BlitSaw s1,s2;
        DspFloatType _out = 0;
        DspFloatType _duty = 0.5;
        DspFloatType sampleRate=44100.0;

        BlitSquare(DspFloatType freq=220,DspFloatType sampleRate=44100) : OscillatorProcessor()
        {
            this->sampleRate = sampleRate;
            block.setFc(10.0f/sampleRate);
            _out = 0;
            _duty = 0.5;
            setFrequency(freq);
            s1.setGain(1);
            s2.setGain(1);
        }
        void setFrequency(DspFloatType f)
        {
            s1.setFrequency(f);
            s2.setFrequency(f);        
        }
        void setDuty(DspFloatType d)
        {
            _duty = d;
        }
        void reset() {
            s1.reset();
            s2.reset();
        }
        enum {
            PORT_FREQ,
            PORT_DUTY,
            PORT_RESET,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_DUTY: setDuty(v); break;
                case PORT_RESET: reset(); break;
                default: printf("No port %d\n", port);
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0) {
            DspFloatType r1 = s1.Tick();        
            s2.setPhaseOffset(s1.getPhase() + _duty*M_PI);
            DspFloatType r2 = s2.Tick();                
            _out = r2-r1;
            DspFloatType x = _out;
            x -= block.process(x);
            return x;
        }

        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
                DspFloatType r1 = s1.Tick();        
                s2.setPhaseOffset(s1.getPhase() + _duty*M_PI);
                DspFloatType r2 = s2.Tick();                
                _out = (r2-r1);
                DspFloatType x = _out;
                x -= block.process(x);
                // this will make it perfectly square but it will most likely alias again
                DspFloatType gibbs;
                if(x < 0) gibbs = - 1.0 + fabs(x);
                else gibbs = 1.0 - x;
                out[i] = x+gibbs;
                if(in) out[i] *= in[i];
            }
        }

        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
            
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,nullptr,input);
        }
    };
}
