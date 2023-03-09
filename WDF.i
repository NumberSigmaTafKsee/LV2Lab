%module WDF
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>


%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

%include "Analog/VAWDFCompressor.hpp"
%include "Analog/VAWDFDiodeClipper.hpp"
%include "Analog/VAWDFPassiveLPF.hpp"
%include "Analog/VAWDFSallenKey.hpp"
