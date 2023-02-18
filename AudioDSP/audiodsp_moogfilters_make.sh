kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_moogfilters.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude  \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o MoogFilters.so $kit/audiodsp_moogfilters_wrap.cxx \
-lstdc++ -lm -lluajit
