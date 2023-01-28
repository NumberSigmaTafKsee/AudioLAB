swig -lua -c++ audiofft.i
gcc -DAUDIOFFT_FFTW3 -O2 -march=native -fPIC -mavx2 -shared -o audiofft.so audiofft_wrap.cxx AudioFFT.cpp -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
