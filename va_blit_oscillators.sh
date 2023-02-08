swig -lua -c++ -Iinclude va_blit_oscillators.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o va_blit_oscillators.so \
va_blit_oscillators_wrap.cxx -lstdc++ -lm -lluajit
