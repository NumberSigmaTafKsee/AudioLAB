swig -lua -c++ -Itcc tcc.i
gcc -Itcc -O2 -fPIC -march=native -mavx2 -shared -o tcc.so tcc_wrap.cxx lib/libtcc.a lib/libtcc1.a -lstdc++ -lm -lluajit -pthread
