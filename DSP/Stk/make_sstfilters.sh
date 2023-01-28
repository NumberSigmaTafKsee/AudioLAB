swig -lua -c++ -Iinclude sstfilters.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -o sstfilters.so sstfilters_wrap.cxx  -lstdc++ -lm -lluajit
