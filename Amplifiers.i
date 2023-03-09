%module Amplifiers
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "ClipFunctions.hpp"
#include "FX/Amplifier.hpp"
#include "FX/Amplifiers.hpp"
#include "FX/Diode.hpp"
#include "FX/DiodeClipper.hpp"
#include "FX/ClipperCircuit.hpp"
#include "FX/ClippingFunctions.hpp"
#include "FX/DistortionFunctions.hpp"
#include "FX/Waveshapers.hpp"
#include "FX/Waveshaping.hpp"
#include "FX/SstWaveshaper.hpp"

%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
//%include "lua_fnptr.i"

%include "SoundObject.hpp"

%include "ClipFunctions.hpp"
%include "FX/Amplifier.hpp"
%include "FX/Amplifiers.hpp"
%include "FX/Diode.hpp"
%include "FX/DiodeClipper.hpp"
%include "FX/ClipperCircuit.hpp"
%include "FX/ClippingFunctions.hpp"
%include "FX/DistortionFunctions.hpp"
%include "FX/Waveshapers.hpp"
%include "FX/Waveshaping.hpp"
%include "FX/SstWaveshaper.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;
