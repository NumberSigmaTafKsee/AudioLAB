kit="Kits/AudioDSP"
kit="."
swig -lua -c++ -Iinclude $kit/audiodsp_fxdsp.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude  \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o FxDSP.so $kit/audiodsp_fxdsp_wrap.cxx \
-lstdc++ -lm -lluajit
