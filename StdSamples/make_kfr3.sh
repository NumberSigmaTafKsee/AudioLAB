swig -lua -c++ -IDSP/Kfr3 kfr3.i
g++ -Iinclude -Isrc -I/usr/local/include/kfr -I../include/Kfr -std=c++17 -fmax-errors=1 -O2 -march=native -fPIC -shared -o kfr3.so kfr3_wrap.cxx -lstdc++ -lluajit -lkfr_dft -lkfr_io -lkfr_capi -lfftw3 -lfftw3f 
