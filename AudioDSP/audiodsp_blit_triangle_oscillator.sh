swig -lua -c++ -Iinclude audiodsp_blit_triangle_oscillator.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o audiodsp_blit_triangle_oscillator.so \
audiodsp_blit_triangle_oscillator.cxx -lstdc++ -lm -lluajit