%module audio_simple_resampler
%{
#define DSPFLOATDOUBLE    
#include "SoundObject.hpp"
#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "audio_simple_resampler.hpp"
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"

#define DSPFLOATDOUBLE    
%include "SoundObject.hpp"

%include "audio_simple_resampler.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

