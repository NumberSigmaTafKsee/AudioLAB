swig -lua -c++ convolution_fftw.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o convolution_fftw.so convolution_fftw_wrap.cxx -lstdc++ -lm -lluajit -lfftw3
