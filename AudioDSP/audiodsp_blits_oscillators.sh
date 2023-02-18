swig -lua -c++ -Iinclude audiodsp_blits_oscillators.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o audiodsp_blits_oscillators.so \
audiodsp_blits_oscillators_wrap.cxx -lstdc++ -lm -lluajit
