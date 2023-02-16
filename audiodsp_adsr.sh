swig -lua -c++ -Iinclude audio_adsr.i
gcc -fmax-errors=1 -fopenmp -pthread -std=c++17 -I. -Iinclude -I.-O2 -fPIC -march=native -mavx2 -shared -fopenmp -pthread -o audio_adsr.so \
audio_adsr_wrap.cxx -lstdc++ -lm -lluajit
