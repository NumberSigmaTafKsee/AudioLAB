swig -lua -c++ -Iinclude/uKfr -Iinclude -Isrc -IKfrDSP kfr2f.i
g++ -Iinclude -Isrc -Iinclude/uKfr -IKfrDSP -std=c++17 -fmax-errors=1 -O2 -march=native -fPIC -shared -o bin/kfr2f.so kfr2f_wrap.cxx lib/libaudiofft.a  -lstdc++ -lluajit -lkfr_dft -lkfr_io -lkfr_capi -lfftw3 -lfftw3f 

