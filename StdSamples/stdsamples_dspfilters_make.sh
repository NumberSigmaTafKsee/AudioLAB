kit="Kits/AudioDSP/"
kit="."
swig -lua -c++ -Iinclude/DSPFilters $kit/audiodsp_lua_dspfilters.i
gcc -fmax-errors=1 -std=c++17 -Iinclude -Iinclude/DSPFilters -O2 -fPIC \
-march=native -mavx2 -shared -o DspFilters.so $kit/audiodsp_DspFilters_wrap.cxx -lstdc++ \
-lm -lluajit -lDSPFilters
