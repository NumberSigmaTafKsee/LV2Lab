swig -lua -c++ -I../AudioLAB/include Amplifiers.i
gcc -fmax-errors=1 -std=c++17 -I. -I../AudioLAB/include -O2 -fPIC -mavx2 -mfma -march=native -shared -o Amplifiers.so Amplifiers_wrap.cxx -lstdc++ -lm -lluajit
