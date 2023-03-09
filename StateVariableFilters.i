%module StateVariableFilters
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Analog/VAAnalogSVF.hpp"

#include "Analog/VAAnalogSVF.hpp"
#include "Analog/VAStateVariableFilter.hpp"
#include "Analog/VAStateVariableFilter1.hpp"
#include "Analog/VAStateVariableFilter2.hpp"
#include "Analog/VAStateVariableFilters.hpp"
#include "Analog/VASVF.hpp"
#include "Analog/VASVFChamberlinFilter.hpp"
#include "Analog/VASVFFilter.hpp"
#include "Analog/VASVFSmoother.hpp"
#include "Analog/VASVSmoothFilter.hpp"
#include "Analog/VASVStateVariableFilter.hpp"
#include "Analog/VAGenSVF.hpp"

%}


%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"

%include "SoundObject.hpp"


// SVF
%include "Analog/VAAnalogSVF.hpp"
%include "Analog/VAGenSVF.hpp"
%rename Analog::Filters::StateVariableFilter::SVFLowpass            VASVFLowPassFilter;
%rename Analog::Filters::StateVariableFilter::SVFBandpass           VASVFBandPassFilter;
%rename Analog::Filters::StateVariableFilter::SVFHighpass           VASVFHighPassFilter;
%rename Analog::Filters::StateVariableFilter::SVFUnitGainBandpass   VASVFUnitGainBandpassFilter;
%rename Analog::Filters::StateVariableFilter::SVFBandShelving       VASVFBandShelvingFilter;
%rename Analog::Filters::StateVariableFilter::SVFNotch              VASVFNotchFilter;
%rename Analog::Filters::StateVariableFilter::SVFPeak               VASVFPeakFilter;
%rename Analog::Filters::StateVariableFilter::SVFAllpass            VASVFAllpassFilter;
%ignore Analog::Filters::StateVariableFilter::resonanceToQ;
%include "Analog/VAStateVariableFilter.hpp"
%include "Analog/VAStateVariableFilter1.hpp"
%rename Analog::Filters::StateVariableFilter2::StateVariableFilter VAStateVariableFilter2;
%include "Analog/VAStateVariableFilter2.hpp"
%rename Analog::Filters::SVF::AnalogSVF VASVF;
%include "Analog/VASVF.hpp"
%include "Analog/VASVFChamberlinFilter.hpp"
%rename Analog::Filters::SVF::StateVariableFilter VASVFFilter;
%include "Analog/VASVFFilter.hpp"
%include "Analog/VASVFSmoother.hpp"
%include "Analog/VASVSmoothFilter.hpp"
