swig -octave -c++ -Iinclude va_blit_square_oscillator.i
mkoctfile -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -mavx2 -mfma -march=native \
-o va_blit_square_oscillator va_blit_square_oscillator_wrap.cxx -lstdc++ -lm -lluajit
