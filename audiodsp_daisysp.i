%module DaisySP
%{

#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

//#include "Daisy.cpp
//#include "Daisy.hpp"
#include "Synthesizer/DaisySP.hpp"
#include "Synthesizer/DaisySPADEnv.hpp"
#include "Synthesizer/DaisySPADSR.hpp"
#include "Synthesizer/DaisySPAllPass.hpp"
#include "Synthesizer/DaisySPAnalogBassDrum.hpp"
#include "Synthesizer/DaisySPAnalogSnareDrum.hpp"
#include "Synthesizer/DaisySPATone.hpp"
#include "Synthesizer/DaisySPAutoWah.hpp"
#include "Synthesizer/DaisySPBalance.hpp"
#include "Synthesizer/DaisySPBiquad.hpp"
#include "Synthesizer/DaisySPBitCrush.hpp"
#include "Synthesizer/DaisySPBLOsc.hpp"
#include "Synthesizer/DaisySPChorus.hpp"
#include "Synthesizer/DaisySPClockedNoise.hpp"
#include "Synthesizer/DaisySPCombFilter.hpp"
#include "Synthesizer/DaisySPCompressor.hpp"
#include "Synthesizer/DaisySPControl.hpp"
#include "Synthesizer/DaisySPCrossFade.hpp"
#include "Synthesizer/DaisySPDCBlock.hpp"
#include "Synthesizer/DaisySPDecimator.hpp"
#include "Synthesizer/DaisySPDelayLine.hpp"
#include "Synthesizer/DaisySPDrip.hpp"
#include "Synthesizer/DaisySPDrums.hpp"
#include "Synthesizer/DaisySPDust.hpp"
#include "Synthesizer/DaisySPDynamics.hpp"
#include "Synthesizer/DaisySPFilters.hpp"
#include "Synthesizer/DaisySPFIR.hpp"
#include "Synthesizer/DaisySPFlanger.hpp"
#include "Synthesizer/DaisySPFM2.hpp"
#include "Synthesizer/DaisySPFold.hpp"
#include "Synthesizer/DaisySPFormantOsc.hpp"
#include "Synthesizer/DaisySPFractalNoise.hpp"
#include "Synthesizer/DaisySPFX.hpp"
#include "Synthesizer/DaisySPGrainlet.hpp"
#include "Synthesizer/DaisySPHarmonicOsc.hpp"
#include "Synthesizer/DaisySPHiHat.hpp"
#include "Synthesizer/DaisySPJFilter.hpp"
#include "Synthesizer/DaisySPJitter.hpp"
#include "Synthesizer/DaisySPKarplusStrong.hpp"
#include "Synthesizer/DaisySPLimiter.hpp"
#include "Synthesizer/DaisySPLine.hpp"
#include "Synthesizer/DaisySPLooper.hpp"
#include "Synthesizer/DaisySPMayTrig.hpp"
#include "Synthesizer/DaisySPMetronome.hpp"
#include "Synthesizer/DaisySPModalVoice.hpp"
#include "Synthesizer/DaisySPMode.hpp"
#include "Synthesizer/DaisySPMoogLadder.hpp"
//#include "Synthesizer/DaisySPNFilt.hpp"
#include "Synthesizer/DaisySPNoise.hpp"
#include "Synthesizer/DaisySPOscillator.hpp"
#include "Synthesizer/DaisySPOscillatorBank.hpp"
#include "Synthesizer/DaisySPOverdrive.hpp"
#include "Synthesizer/DaisySPParticle.hpp"
#include "Synthesizer/DaisySPPhaser.hpp"
#include "Synthesizer/DaisySPPhasor.hpp"
//#include "Synthesizer/DaisySPPhysicalModels.hpp"
#include "Synthesizer/DaisySPPitchShifter.hpp"
#include "Synthesizer/DaisySPPluck.hpp"
#include "Synthesizer/DaisySPPolyPluck.hpp"
#include "Synthesizer/DaisySPPort.hpp"
//#include "Synthesizer/DaisySPResonator.hpp"
#include "Synthesizer/DaisySPReverbSC.hpp"
#include "Synthesizer/DaisySPSampleHold.hpp"
#include "Synthesizer/DaisySPSampleRateReducer.hpp"
#include "Synthesizer/DaisySPSmoothRandom.hpp"
#include "Synthesizer/DaisySPStringVoice.hpp"
#include "Synthesizer/DaisySPSVF.hpp"
#include "Synthesizer/DaisySPSynthBassDrum.hpp"
#include "Synthesizer/DaisySPSynthesis.hpp"
#include "Synthesizer/DaisySPSynthSnareDrum.hpp"
#include "Synthesizer/DaisySPTone.hpp"
#include "Synthesizer/DaisySPTremolo.hpp"
#include "Synthesizer/DaisySPUtils.hpp"
#include "Synthesizer/DaisySPVariableSaw.hpp"
#include "Synthesizer/DaisySPVariableShapeOsc.hpp"
#include "Synthesizer/DaisySPVOSIM.hpp"
#include "Synthesizer/DaisySPWaveFolder.hpp"
#include "Synthesizer/DaisySPWhiteNoise.hpp"
#include "Synthesizer/DaisySPZOscillator.hpp"
%}


%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "Synthesizer/DaisySP.hpp"
%include "Synthesizer/DaisySPADEnv.hpp"
%include "Synthesizer/DaisySPADSR.hpp"
%include "Synthesizer/DaisySPAllPass.hpp"
%include "Synthesizer/DaisySPAnalogBassDrum.hpp"
%include "Synthesizer/DaisySPAnalogSnareDrum.hpp"
%include "Synthesizer/DaisySPATone.hpp"
%include "Synthesizer/DaisySPAutoWah.hpp"
%include "Synthesizer/DaisySPBalance.hpp"
%include "Synthesizer/DaisySPBiquad.hpp"
%include "Synthesizer/DaisySPBitCrush.hpp"
%include "Synthesizer/DaisySPBLOsc.hpp"
%include "Synthesizer/DaisySPChorus.hpp"
%include "Synthesizer/DaisySPClockedNoise.hpp"
%include "Synthesizer/DaisySPCombFilter.hpp"
%include "Synthesizer/DaisySPCompressor.hpp"
//%include "Synthesizer/DaisySPControl.hpp"
%include "Synthesizer/DaisySPCrossFade.hpp"
%include "Synthesizer/DaisySPDCBlock.hpp"
%include "Synthesizer/DaisySPDecimator.hpp"
%include "Synthesizer/DaisySPDelayLine.hpp"
%include "Synthesizer/DaisySPDrip.hpp"
//%include "Synthesizer/DaisySPDrums.hpp"
%include "Synthesizer/DaisySPDust.hpp"
//%include "Synthesizer/DaisySPDynamics.hpp"
//%include "Synthesizer/DaisySPFilters.hpp"
%include "Synthesizer/DaisySPFIR.hpp"
%include "Synthesizer/DaisySPFlanger.hpp"
%include "Synthesizer/DaisySPFM2.hpp"
%include "Synthesizer/DaisySPFold.hpp"
%include "Synthesizer/DaisySPFormantOsc.hpp"
%include "Synthesizer/DaisySPFractalNoise.hpp"
//%include "Synthesizer/DaisySPFX.hpp"
%include "Synthesizer/DaisySPGrainlet.hpp"
%include "Synthesizer/DaisySPHarmonicOsc.hpp"
%include "Synthesizer/DaisySPHiHat.hpp"
%include "Synthesizer/DaisySPJFilter.hpp"
%include "Synthesizer/DaisySPJitter.hpp"
%include "Synthesizer/DaisySPKarplusStrong.hpp"
%include "Synthesizer/DaisySPLimiter.hpp"
%include "Synthesizer/DaisySPLine.hpp"
%include "Synthesizer/DaisySPLooper.hpp"
%include "Synthesizer/DaisySPMayTrig.hpp"
%include "Synthesizer/DaisySPMetronome.hpp"
%include "Synthesizer/DaisySPModalVoice.hpp"
%include "Synthesizer/DaisySPMode.hpp"
%include "Synthesizer/DaisySPMoogLadder.hpp"
//%include "Synthesizer/DaisySPNFilt.hpp"
//%include "Synthesizer/DaisySPNoise.hpp"
%include "Synthesizer/DaisySPOscillator.hpp"
%include "Synthesizer/DaisySPOscillatorBank.hpp"
%include "Synthesizer/DaisySPOverdrive.hpp"
%include "Synthesizer/DaisySPParticle.hpp"
%include "Synthesizer/DaisySPPhaser.hpp"
%include "Synthesizer/DaisySPPhasor.hpp"
//%include "Synthesizer/DaisySPPhysicalModels.hpp"
%include "Synthesizer/DaisySPPitchShifter.hpp"
%include "Synthesizer/DaisySPPluck.hpp"
%include "Synthesizer/DaisySPPolyPluck.hpp"
%include "Synthesizer/DaisySPPort.hpp"
//%include "Synthesizer/DaisySPResonator.hpp"
%include "Synthesizer/DaisySPReverbSC.hpp"
%include "Synthesizer/DaisySPSampleHold.hpp"
%include "Synthesizer/DaisySPSampleRateReducer.hpp"
%include "Synthesizer/DaisySPSmoothRandom.hpp"
%include "Synthesizer/DaisySPStringVoice.hpp"
%include "Synthesizer/DaisySPSVF.hpp"
%include "Synthesizer/DaisySPSynthBassDrum.hpp"
//%include "Synthesizer/DaisySPSynthesis.hpp"
%include "Synthesizer/DaisySPSynthSnareDrum.hpp"
%include "Synthesizer/DaisySPTone.hpp"
%include "Synthesizer/DaisySPTremolo.hpp"
//%include "Synthesizer/DaisySPUtils.hpp"
%include "Synthesizer/DaisySPVariableSaw.hpp"
%include "Synthesizer/DaisySPVariableShapeOsc.hpp"
%include "Synthesizer/DaisySPVOSIM.hpp"
%include "Synthesizer/DaisySPWaveFolder.hpp"
%include "Synthesizer/DaisySPWhiteNoise.hpp"
%include "Synthesizer/DaisySPZOscillator.hpp"


%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%inline %{
    const int BufferSize = 256;
    Default noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
%}
