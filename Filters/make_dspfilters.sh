swig -lua -c++ -Iinclude -Iinclude/Filters/DSPFilters audiodsp_dspfilters.i
gcc -fmax-errors=1 -std=c++17 -Iinclude/Filters/DSPFilters -Iinclude -O2 -fPIC \
-march=native -mavx2 -shared -o dspfilters.so audiodsp_dspfilters_wrap.cxx -lstdc++ -lm -lluajit \
-lDSPFilters
