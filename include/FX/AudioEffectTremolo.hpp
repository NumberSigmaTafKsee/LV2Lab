#pragma once

namespace FX
{
    ///////////////////////////////////////////////////////////////////////////////
    // Tremolo
    ///////////////////////////////////////////////////////////////////////////////
    struct Tremolo : public StereoFXProcessor
    {
        Waveform Tremolo_LFO[2];
        DspFloatType    Tremolo_ModDepth[2];
        DspFloatType    Tremolo_LfoPhase[2];
        DspFloatType    Tremolo_LfoFreqHz[2];
        DspFloatType    sampleRate,inverseSampleRate;
        int      Tremolo_NumChannels;
        
        Tremolo(DspFloatType sr=44100,Waveform lfo=kWaveformSine, DspFloatType modDepth=0.95, DspFloatType freq=2.1, size_t nChannels=2)
        : StereoFXProcessor()
        {
            sampleRate = sr;
            inverseSampleRate = 1.0/sr;
            setData<DspFloatType>(Tremolo_ModDepth,modDepth);
            setData<DspFloatType>(Tremolo_LfoFreqHz,freq);
            setData<Waveform>(Tremolo_LFO,lfo);
            setData<DspFloatType>(Tremolo_LfoPhase,0);
            Tremolo_NumChannels = nChannels;
        }
        enum {
            PORT_LFO1_WAVEFORM,
            PORT_LFO2_WAVEFORM,
            PORT_MODDEPTH1,
            PORT_MODDEPTH2,
            PORT_LFOPHASE1,
            PORT_LFOPHASE2,
            PORT_LFOFREQ1,
            PORT_LFOFREQ2,            
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_LFO1_WAVEFORM: Tremolo_LFO[0] = (Waveform)v; break;
                case PORT_LFO2_WAVEFORM: Tremolo_LFO[1] = (Waveform)v; break;
                case PORT_MODDEPTH1: Tremolo_ModDepth[0] = v; break;
                case PORT_MODDEPTH2: Tremolo_ModDepth[1] = v; break;
                case PORT_LFOPHASE1: Tremolo_LfoPhase[0] = v; break;
                case PORT_LFOPHASE2: Tremolo_LfoPhase[1] = v; break;
                case PORT_LFOFREQ1: Tremolo_LfoFreqHz[0] = v; break;
                case PORT_LFOFREQ2: Tremolo_LfoFreqHz[1] = v; break;
                default: printf("No port %d\n",port);
            }
        }
        // Process one buffer ("block") of data
        void ProcessBlock (size_t n, DspFloatType ** inputs, DspFloatType ** outputs)
        {
            Undenormal denormal;

            // apply the same modulation to all input channels for which there is an output channel
            
            for (int channelIndex = 0; channelIndex < Tremolo_NumChannels; channelIndex++)
            {
                DspFloatType phi = Tremolo_LfoPhase[channelIndex];
                DspFloatType deltaPhi = DspFloatType(Tremolo_LfoFreqHz[channelIndex] * inverseSampleRate);

                // restart the phase sequence
                phi = Tremolo_LfoPhase[channelIndex];

                DspFloatType* pIn  = inputs[channelIndex]  ;
				DspFloatType* pOut = outputs[channelIndex] ;
				#pragma omp simd aligned(pIn,pOut,inputs,outputs)
                for (int i = 0; i < n; i++)
                {
                    DspFloatType modAmount = LFO_GetSample(phi, Tremolo_LFO[channelIndex]);
                    DspFloatType x = (1.0f - Tremolo_ModDepth[channelIndex] * modAmount);
                    
                    *pOut++ = *pIn++ * x;                
                    // Update LFO phase, keeping in range [0, 1]
                    phi += deltaPhi;
                    while (phi >= 1.0) phi -= 1.0;
                }
                // update the main LFO phase state variable, ready for the next processBlock() call
                Tremolo_LfoPhase[channelIndex] = phi;    
            }        
        }
        void ProcessInplace(size_t n, DspFloatType ** buffer)
        {
            ProcessBlock(n,buffer,buffer);
        }
    };
}
