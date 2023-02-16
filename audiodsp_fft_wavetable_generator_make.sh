swig -lua -c++ -Iinclude audio_fft_wavetable_generator.i
gcc -fmax-errors=1 -std=c++17 -I. -Iinclude -I. -O2 -fPIC -march=native -mavx2 -shared -o audio_fft_wavetable_generator.so \
audio_fft_wavetable_generator_wrap.cxx -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
