swig -lua -c++ -Iinclude audio_wavetable.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I. -O2 -fPIC -march=native -mavx2 -shared -o audio_wavetable.so \
audio_wavetable_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
