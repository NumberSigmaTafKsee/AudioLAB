swig -lua -c++ -Iinclude/Filters/DSPFilters audiodsp_DspFilters.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -Iinclude/Filters/DSPFilters -O2 -fPIC \
-march=native -mavx2 -shared -o DspFilters.so audiodsp_DspFilters_wrap.cxx -lstdc++ \
-lm -lluajit -lDSPFilters
