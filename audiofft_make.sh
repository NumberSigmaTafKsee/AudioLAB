swig -lua -c++ -Iinclude/AudioFFT Kits/AudioDSP/audiodsp_audiofft.i
gcc -Iinclude/AudioFFT -DAUDIOFFT_FFTW3 -O2 -march=native -fPIC -mavx2 -shared -o audiofft.so \
Kits/AudioDSP/audiodsp_audiofft_wrap.cxx include/AudioFFT/AudioFFT.cpp -lstdc++ -lm -lluajit -lfftw3 -lfftw3f
