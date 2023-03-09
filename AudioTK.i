%module AudioTK
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "AudioTK/ATK.hpp"
#include "AudioTK/ATKAdaptiveFilters.hpp"
#include "AudioTK/ATKAttackRelease.hpp"
#include "AudioTK/ATKAttackReleaseHysterisis.hpp"
#include "AudioTK/ATKBesselFilters.hpp"
#include "AudioTK/ATKBlockLMSFilter.hpp"
#include "AudioTK/ATKButterworthFilters.hpp"
#include "AudioTK/ATKChamberlinFilter.hpp"
#include "AudioTK/ATKChebyshev1Filter.hpp"
#include "AudioTK/ATKChebyshev2Filter.hpp"
#include "AudioTK/ATKDelays.hpp"
#include "AudioTK/ATKDistortionProcessors.hpp"
#include "AudioTK/ATKDynamicProcessors.hpp"
#include "AudioTK/ATKEqProcessors.hpp"
#include "AudioTK/ATKExpander.hpp"
#include "AudioTK/ATKFeedbackDelayNetwork.hpp"
#include "AudioTK/ATKFIRFilter.hpp"
#include "AudioTK/ATKFixedDelayLine.hpp"
#include "AudioTK/ATKFollowerTransistorClassAProcessor.hpp"
#include "AudioTK/ATKGainColoredCompressor.hpp"
#include "AudioTK/ATKGainColoredExpander.hpp"
#include "AudioTK/ATKGainCompressor.hpp"
//#include "AudioTK/ATKGainFilter.hpp"
#include "AudioTK/ATKGainLimiter.hpp"
#include "AudioTK/ATKGainMaxColoredExpander.hpp"
#include "AudioTK/ATKGainSwell.hpp"
#include "AudioTK/ATKIIRFilter.hpp"
#include "AudioTK/ATKLinkwitzReillyFilters.hpp"
#include "AudioTK/ATKLMSFilter.hpp"
#include "AudioTK/ATKMaxCompressor.hpp"
#include "AudioTK/ATKMaxExpander.hpp"
#include "AudioTK/ATKMultipleUniversalFixedDelayLine.hpp"

#include "AudioTK/ATKPowerFilter.hpp"
//#include "AudioTK/ATKPreampProcessors.hpp"
#include "AudioTK/ATKRBJFilters.hpp"
#include "AudioTK/ATKRelativePowerFilter.hpp"
#include "AudioTK/ATKRemezFilter.hpp"

#include "AudioTK/ATKReverbProcessors.hpp"
#include "AudioTK/ATKRIAAFilters.hpp"
#include "AudioTK/ATKRLSFilter.hpp"
#include "AudioTK/ATKSecondOrderFilters.hpp"
#include "AudioTK/ATKSecondOrderSVFFilters.hpp"
#include "AudioTK/ATKTimeVaryingFilters.hpp"
#include "AudioTK/ATKTimeVaryingSVFFilters.hpp"
#include "AudioTK/ATKToneFilters.hpp"
#include "AudioTK/ATKToneStackFilter.hpp"
#include "AudioTK/ATKToolProcessors.hpp"
#include "AudioTK/ATKTransistorClassAProcessor.hpp"
#include "AudioTK/ATKTriode2Processor.hpp"
#include "AudioTK/ATKTriodeProcessor.hpp"
#include "AudioTK/ATKUniversalFixedDelayLine.hpp"
#include "AudioTK/ATKUniversalVariableDelayLine.hpp"
#include "AudioTK/ATKVariableDelayLine.hpp"


%}
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

%include "AudioTK/ATK.hpp"
%include "AudioTK/ATKAttackRelease.hpp"
%include "AudioTK/ATKAttackReleaseHysterisis.hpp"
%include "AudioTK/ATKBesselFilters.hpp"
%include "AudioTK/ATKBlockLMSFilter.hpp"
%include "AudioTK/ATKButterworthFilters.hpp"
%include "AudioTK/ATKChamberlinFilter.hpp"
%include "AudioTK/ATKChebyshev1Filter.hpp"
%include "AudioTK/ATKChebyshev2Filter.hpp"
%include "AudioTK/ATKExpander.hpp"
%include "AudioTK/ATKFeedbackDelayNetwork.hpp"
%include "AudioTK/ATKFIRFilter.hpp"
%include "AudioTK/ATKFixedDelayLine.hpp"
%include "AudioTK/ATKFollowerTransistorClassAProcessor.hpp"
%include "AudioTK/ATKGainColoredCompressor.hpp"
%include "AudioTK/ATKGainColoredExpander.hpp"
%include "AudioTK/ATKGainCompressor.hpp"
//%include "AudioTK/ATKGainFilter.hpp"
%include "AudioTK/ATKGainLimiter.hpp"
%include "AudioTK/ATKGainMaxColoredExpander.hpp"
%include "AudioTK/ATKGainSwell.hpp"
%include "AudioTK/ATKIIRFilter.hpp"
%include "AudioTK/ATKLinkwitzReillyFilters.hpp"
%include "AudioTK/ATKLMSFilter.hpp"
%include "AudioTK/ATKMaxCompressor.hpp"
%include "AudioTK/ATKMaxExpander.hpp"
%include "AudioTK/ATKMultipleUniversalFixedDelayLine.hpp"
%include "AudioTK/ATKPowerFilter.hpp"
%include "AudioTK/ATKPreampProcessors.hpp"
%include "AudioTK/ATKRBJFilters.hpp"
%include "AudioTK/ATKRelativePowerFilter.hpp"
//%include "AudioTK/ATKRemezFilter.hpp"
%include "AudioTK/ATKReverbProcessors.hpp"
%include "AudioTK/ATKRIAAFilters.hpp"
%include "AudioTK/ATKRLSFilter.hpp"
%include "AudioTK/ATKSecondOrderFilters.hpp"
%include "AudioTK/ATKSecondOrderSVFFilters.hpp"
%include "AudioTK/ATKTimeVaryingFilters.hpp"
%include "AudioTK/ATKTimeVaryingSVFFilters.hpp"
%include "AudioTK/ATKToneFilters.hpp"
%include "AudioTK/ATKToneStackFilter.hpp"
%include "AudioTK/ATKToolProcessors.hpp"
%include "AudioTK/ATKTransistorClassAProcessor.hpp"
%include "AudioTK/ATKTriode2Processor.hpp"
%include "AudioTK/ATKTriodeProcessor.hpp"
%include "AudioTK/ATKUniversalFixedDelayLine.hpp"
%include "AudioTK/ATKUniversalVariableDelayLine.hpp"
%include "AudioTK/ATKVariableDelayLine.hpp"


%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;


