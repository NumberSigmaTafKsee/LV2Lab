swig -lua -c++ audiodsp_kissfft.i
gcc -I/usr/local/include/kissfft -L/usr/local/lib -O2 -march=native -mavx2 -fPIC \
-shared -o kissfft.so audiodsp_kissfft_wrap.cxx -lstdc++ -lm -lluajit-5.1 -lkissfft-float
