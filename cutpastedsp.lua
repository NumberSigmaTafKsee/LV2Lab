template = [[
~module %s
~{
typedef float DspFloatType;
#include "SoundObject.hpp"
#include "%s"
~}
typedef float DspFloatType;
~include "stdint.i"
~include "std_math.i"
~include "std_vector.i"
~include "std_list.i"
~include "std_map.i"
~include "lua_fnptr.i"

typedef float DspFloatType;
~include "SoundObject.hpp"
~include "%s"

~template(float_vector) std::vector<float>;
~template(double_vector) std::vector<double>;

~template(complex_float_vector) std::vector<std::complex<float>>;
~template(complex_double_vector) std::vector<std::complex<double>>;

]]

make = [[
    swig -lua -c++ -Iinclude %s.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o %s.so %s_wrap.cxx -lstdc++ -lm -lluajit 
    ]]
    
function make_files(dir,hdr)
    x = string.format(template,hdr, hdr .. ".hpp", hdr .. ".hpp")
    x = x:gsub("~","%%")    
    f = io.open(hdr .. ".swg",'w')
    f:write(x)
    f:close()

    x = string.format(make,hdr,hdr,hdr)
    f = io.open("make_" .. hdr .. ".sh",'w')
    f:write(x)
    f:close()
    os.execute("sh make_" .. hdr .. ".sh")
end


for i=2,#arg  do	
	file=arg[i]    
	if(file:find("hpp") ~= 0) then		
		make_files(arg[1],file:gsub(".hpp",""))    
	end
end
