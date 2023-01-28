swig -lua -c++ -Iinclude src/fir1.i
gcc -Iinclude -O2 -march=native -mavx2 -fPIC -shared -o fir1.so src/fir1_wrap.cxx -lstdc++ -lm -lluajit -Lbin -lfir
