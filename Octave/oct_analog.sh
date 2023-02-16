swig -octave -c++ -Iinclude Analog.i
mkoctfile -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -mavx2 -mfma -march=native \
-o Analog Analog_wrap.cxx -lstdc++ -lm -lluajit
