swig -octave -c++ -Iinclude va_blit_oscillators.i
mkoctfile -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -mavx2 -mfma -march=native \
-o va_blit_oscillators va_blit_oscillators_wrap.cxx -lstdc++ -lm -lluajit
