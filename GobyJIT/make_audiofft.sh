swig -lua -c++ src/audiofft.i
gcc -DAUDIOFFT_FFTW3 -O2 -march=native -fPIC -mavx2 -shared -o audiofft.so src/audiofft_wrap.cxx src/AudioFFT.cpp -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
