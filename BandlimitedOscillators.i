%module BandlimitedOscillators
%{

#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Analog/VABlitSawOscillator.hpp"
#include "Analog/VABlitSquareOscillator.hpp"
#include "Analog/VABlitTriangleOscillator.hpp"
#include "Analog/VablitOscillators.hpp"
#include "Analog/VADPWSawOscillator.hpp"    
#include "Analog/VADPWPulseOscillator.hpp"
#include "Analog/VADPWTriangleOscillator.hpp"
#include "Analog/VAMinBLEP.hpp"
#include "Analog/VAPolyBLEPOscillator.hpp"

%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"

%include "SoundObject.hpp"

%include "Analog/VABlitSawOscillator.hpp"
%include "Analog/VABlitSquareOscillator.hpp"
%include "Analog/VABlitTriangleOscillator.hpp"
%include "Analog/VablitOscillators.hpp"
%include "Analog/VADPWSawOscillator.hpp"    
%include "Analog/VADPWPulseOscillator.hpp"
%include "Analog/VADPWTriangleOscillator.hpp"
%include "Analog/VAMinBLEP.hpp"
%ignore Analog::Oscillators::PolyBLEPOsc::blep;
%ignore Analog::Oscillators::PolyBLEPOsc::blamp;
%include "Analog/VAPolyBLEPOscillator.hpp"
