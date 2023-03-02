swig -octave -c++ -Iinclude Amplifiers.i
mkoctfile -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -mavx2 -mfma -march=native \
-o Amplifier Amplifier_wrap.cxx -lstdc++ -lm -lluajit
