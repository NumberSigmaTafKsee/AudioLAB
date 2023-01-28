swig -lua -c++ DspFilters.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DspFilters.so DspFilters_wrap.cxx -lstdc++ -lm -lluajit -lDSPFilters
