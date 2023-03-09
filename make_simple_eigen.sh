swig -lua -c++ -Iinclude simple_eigen.i
gcc -Wfatal-errors -Iinclude -fPIC -O2 -march=native -mavx2 -mfma -shared \
-o se1.so simple_eigen_wrap.cxx -lstdc++ -lm -lluajit
