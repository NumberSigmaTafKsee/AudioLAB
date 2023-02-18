swig -lua -c++ -Iinclude Kits/AudioDSP/audiodsp_vaoscillators.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -IKits/AudioDSP \
-O2 -fPIC -mavx2 -mfma -march=native -shared \
-o VAOscillators.so Kits/AudioDSP/audiodsp_vaoscillators_wrap.cxx \
-lstdc++ -lm -lluajit
