swig -python -c++ audiofft.i
gcc -DAUDIOFFT_FFTW3 -I/usr/local/include/python3.9 -O2 -march=native -fPIC -mavx2 -shared -o _audiofft.so audiofft_wrap.cxx AudioFFT.cpp -lstdc++ -lm -lpython3.9 -lfftw3 -lfftw3f
