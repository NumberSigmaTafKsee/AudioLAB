swig -lua -c++ dspfilters.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o dspfilters.so dspfilters_wrap.cxx -lstdc++ -lm -lluajit -lDSPFilters
