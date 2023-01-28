swig -lua -c++ src/adsr.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o adsr.so src/adsr_wrap.cxx -lm -lluajit
