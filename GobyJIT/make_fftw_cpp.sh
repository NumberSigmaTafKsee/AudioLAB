swig -lua -c++ fftw2d_cpp.i
swig -lua -c++ fftw_cpp.i
gcc -Ifftw_cpp -O2 -fPIC -march=native -mavx2 -shared -o fftw2d_cpp.so fftw2d_cpp_wrap.cxx -lstdc++ -lm -lluajit -lfftw3
gcc -Ifftw_cpp -O2 -fPIC -march=native -mavx2 -shared -o fftw_cpp.so fftw_cpp_wrap.cxx -lstdc++ -lm -lluajit -lfftw3
