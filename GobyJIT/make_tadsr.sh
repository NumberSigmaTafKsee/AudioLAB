swig -lua -c++ tadsr.i
gcc -Iinclude -fmax-errors=1 -O2 -fPIC -march=native -mavx2 -shared -o TADSR.so tadsr_wrap.cxx -lstdc++ -lm -lluajit
