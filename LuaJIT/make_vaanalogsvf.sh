swig -lua -c++ vaanalogsvf.i
gcc -std=c++17 -I. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o vaanalogsvf.so vaanalogsvf_wrap.cxx -lstdc++ -lm -lluajit
