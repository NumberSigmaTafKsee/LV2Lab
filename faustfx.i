%module faustfx
%{
#include "Faust/FaustFX.hpp"
#include <vector>
%}

%include "std_vector.i"
%include "std_string.i"

%include "Faust/FaustFX.hpp"

%template(float_vector) std::vector<float>;

%inline %{
std::list<GUI*> GUI::fGuiList;
ztimedmap GUI::gTimedZoneMap;
%}
