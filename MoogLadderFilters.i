%module MoogLadderFilters
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Analog/VAMoogLadderFilters.hpp"
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"

%include "SoundObject.hpp"

%include "Analog/VAMoogLadderFiltersBase.hpp"
%include "Analog/VAMoogLadderRBJFilter.hpp"
%include "Analog/VAMoogLadderStilson.hpp"
%include "Analog/VAMoogLadderSimplified.hpp"
%include "Analog/VAMoogLadderRungeKutta.hpp"
%include "Analog/VAMoogLadderOberheim.hpp"
%include "Analog/VAMoogLadderMusicDSP.hpp"
%include "Analog/VAMoogLadderMicrotracker.hpp"
%include "Analog/VAMoogLadderKrajeski.hpp"
%include "Analog/VAMoogLadderImproved.hpp"
%include "Analog/VAMoogLadderHouvilainen.hpp"
%include "Analog/VAMoogLadderFilters.hpp"
