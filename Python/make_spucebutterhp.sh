swig -lua -c++ iirspucebutterhp.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o spucebutterhp.so iirspucebutterhp_wrap.cxx -lstdc++ -lm -lluajit -lspuce
