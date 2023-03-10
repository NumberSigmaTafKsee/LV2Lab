gcc -std=c++17 -I../AudioLAB/include -fPIC -O2 -march=native -mavx2 -mfma -shared -o lv2luajit.so lv2luajit.cpp -lstdc++ -lm -llv2-plugin -llv2-gui -lluajit
mv lv2luajit.so LV2/lv2luajit.lv2
