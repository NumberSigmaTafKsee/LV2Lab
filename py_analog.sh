swig -python -c++ -Iinclude Analog.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I/usr/local/include/python3.9 \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o _Analog.so Analog_wrap.cxx \
-lstdc++ -lm -L/usr/local/lib/python3.9 -lpython3.9 -lutil -ldl -lm -lpthread
