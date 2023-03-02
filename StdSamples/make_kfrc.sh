swig -lua -c++ -I/usr/local/include/kfr -Iinclude/uKfr -Iinclude -Isrc -IKfrDSP kfrc.i
g++ -Iinclude -Isrc -I/usr/local/include/kfr -Iinclude/uKfr -IKfrDSP -std=c++17 -fmax-errors=1 -O2 -march=native -fPIC -shared -o bin/kfrc.so kfrc_wrap.cxx lib/libaudiofft.a  -lstdc++ -lluajit -lkfr_dft -lkfr_io -lkfr_capi -lfftw3 -lfftw3f 

