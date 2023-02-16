swig -octave -c++ -Iinclude OctaveAnalog.i
mkoctfile -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -mavx2 -mfma -march=native \
-o Analog OctaveAnalog_wrap.cxx -lstdc++ -lm -lluajit
