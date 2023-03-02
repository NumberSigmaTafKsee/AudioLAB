swig -lua -c++ -Iinclude/uKfr -Iinclude -Isrc kfr2.i
g++ -Iinclude -Isrc -Iinclude -std=c++17 -fmax-errors=1 -O2 -march=native -fPIC -shared -o kfr2.so kfr2_wrap.cxx lib/libaudiofft.a  -lstdc++ -lluajit -lkfr_dft -lkfr_io -lkfr_capi -lfftw3 -lfftw3f 

