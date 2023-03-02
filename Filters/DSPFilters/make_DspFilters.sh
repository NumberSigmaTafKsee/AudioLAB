swig -lua -c++ -I../../include DspFilters.i
gcc -fmax-errors=1 -std=c++17 -I../../include -O2 -fPIC -march=native -mavx2 -shared -o DspFilters.so DspFilters_wrap.cxx -lstdc++ -lm -lluajit -lDSPFilters
