swig -Iinclude -lua -c++ audiodsp_noisey.i
gcc -Iinclude -fmax-errors=1 -O2 -fPIC -shared -o noisey.so audiodsp_noisey_wrap.cxx -lstdc++ -lm -lluajit
