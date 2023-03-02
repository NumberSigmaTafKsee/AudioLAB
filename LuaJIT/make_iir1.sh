swig -lua -c++ -Iinclude src/iir1.i
gcc -Iinclude -O2 -march=native -mavx2 -fPIC -shared -o iir1.so src/iir1_wrap.cxx -lstdc++ -lm -lluajit -Lbin -lfir
