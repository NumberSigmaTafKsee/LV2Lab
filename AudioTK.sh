swig -lua -c++ -I../AudioLAB/include -I../AudioLAB AudioTK.i
gcc -fmax-errors=1 -std=c++17 -I../AudioLAB/include -I../AudioLAB -I/usr/local/include/lilv-0 -O2 -fPIC -mavx2 -mfma -march=native -shared \
-o AudioTK.so AudioTK_wrap.cxx -lstdc++ -lm -lluajit -lATKCore -lATKAdaptive -lATKDelay -lATKDistortion -lATKDynamic -lATKEQ -lATKIO -lATKPreamplifier \
-lATKReverberation -lATKSpecial -lATKTools -lATKUtility
