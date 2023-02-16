swig -lua -c++ -Iinclude vafilters.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o vafilters.so vafilters_wrap.cxx -lstdc++ -lm -lluajit
