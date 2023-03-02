swig -lua -c++ -Isrc -Iinclude src/keyfinder.i
gcc -Isrc -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o keyfinder.so src/keyfinder_wrap.cxx -lstdc++ -lm -lluajit -Lbin -lkeyfinder
