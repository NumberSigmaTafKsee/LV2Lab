%module Diodes
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>


#include "Analog/VADiodeLadderFilter2.hpp"
#include "Analog/VirtualAnalogDiodeLadderFilter.hpp"
#include "Analog/VADiode.hpp"
#include "Analog/VADiodeClipper.hpp"
#include "Analog/VADiodeLadderFilter.cpp
#include "Analog/VADiodeLadderFilter.hpp"
#include "Analog/VADiodeSimulator.hpp"
#include "Analog/VAVCS3DiodeFilter.hpp"
#include "Analog/VAVCS3Filter.hpp"

%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"

%include "SoundObject.hpp"


%include "Analog/VADiodeLadderFilter2.hpp"
%include "Analog/VADiode.hpp"
%include "Analog/VADiodeClipper.hpp"
%include "Analog/VADiodeLadderFilter.cpp
%include "Analog/VADiodeLadderFilter.hpp"
%include "Analog/VADiodeSimulator.hpp"
%include "Analog/VAVCS3DiodeFilter.hpp"
%include "Analog/VAVCS3Filter.hpp"

//%rename Analog::Filters::VirtualAnalogDiodeLadderFilter VADiodeLadderFilter;
//%include "Analog/VirtualAnalogDiodeLadderFilter.hpp"

