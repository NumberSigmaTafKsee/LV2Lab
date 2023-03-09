%module VoltageControlled
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Analog/VCA.hpp"
#include "Analog/VCF.hpp"
#include "Analog/VCO.hpp"

%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"

%include "SoundObject.hpp"

%include "Analog/VCA.hpp"
%include "Analog/VCF.hpp"
%include "Analog/VCO.hpp"

//%include "Analog/VoltageControlledFilter.hpp"
//%include "Analog/VoltageControlledOscillator.hpp"
