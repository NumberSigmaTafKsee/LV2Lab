swig -lua -c++ -I../AudioLAB/include -Iinclude audiodsp_stk.i 
gcc -D__UNIX_JACK__ -D__LINUX_PULSE__ -DRAWWAVE_PATH=../AudioLAB/Data/rawwaves \
-fmax-errors=1 -std=c++17 -I../AudioLAB/include -I../AudioLAB -O2 -fPIC -mavx2 -mfma -march=native -shared -o Stk.so audiodsp_stk_wrap.cxx \
../AudioLAB/lib/libstk.a -lstdc++ -lm -lluajit
