swig -lua -c++ -I/usr/include src/fftw3.i
gcc -fopenmp -O2 -fPIC -march=native -mavx2 -shared -o fftw3.so src/fftw3_wrap.cxx -lstdc++ -lm -lluajit -lfftw3_omp -lfftw3f_omp -lfftw3l_omp -lfftw3q_omp -lfftw3_threads
