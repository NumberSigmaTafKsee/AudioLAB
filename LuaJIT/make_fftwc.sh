swig -lua -c++ fftw_convolve.i
gcc -fpermissive -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o fftw_convolve.so fftw_convolve_wrap.cxx fftw_wrappers.cc -lstdc++ -lm -lluajit -lfftw3f
