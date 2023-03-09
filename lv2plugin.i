%module lv2plugin
%{
#include "SoundObject.hpp"
#include "LV2/LV2Plugin.hpp"
#include <vector>
%}

%include "std_vector.i"
%include "std_string.i"

%include "SoundObject.hpp"
%include "LV2/LV2Plugin.hpp"

%template(float_vector) std::vector<float>;

