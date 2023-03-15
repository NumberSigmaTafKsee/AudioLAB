swig -lua -c++ -I../../include dspfilters.i
gcc -fmax-errors=1 -std=c++17 -I../../include -O2 -fPIC -march=native -mavx2 -shared -o dspfilters.so dspfilters_wrap.cxx -lstdc++ -lm -lluajit -lDSPFilters
