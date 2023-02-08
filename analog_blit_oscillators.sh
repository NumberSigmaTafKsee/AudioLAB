swig -lua -c++ -Iinclude analog_blit_oscillators.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o analog_blit_oscillators.so \
analog_blit_oscillators_wrap.cxx -lstdc++ -lm -lluajit
