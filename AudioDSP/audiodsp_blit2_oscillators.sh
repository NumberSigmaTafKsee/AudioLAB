swig -octave -c++ -Iinclude audiodsp_blit2_oscillators.i
mkoctfile -Wfatal-errors -fmax-errors=1 -std=c++17 -I. -Iinclude -O2 -fPIC -mavx2 -mfma -march=native -fopenmp -pthread \
-o audiodsp_blit2_oscillators.so audiodsp_blit2_oscillators_wrap.cxx -lstdc++ -lm -lluajit
