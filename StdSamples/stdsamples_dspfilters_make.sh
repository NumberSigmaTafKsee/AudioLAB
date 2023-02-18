kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude -Iinclude/DSPFilters $kit/audiodsp_dspfilters.i
gcc -fmax-errors=1 -std=c++17 -Iinclude/DSPFilters -Iinclude -O2 -fPIC \
-march=native -mavx2 -shared -o dspfilters.so $kit/audiodsp_dspfilters_wrap.cxx -lstdc++ -lm -lluajit \
-lDSPFilters
