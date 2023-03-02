swig -lua -c++ iirspucebutterbp.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o spucebutterbp.so iirspucebutterbp_wrap.cxx -lstdc++ -lm -lluajit -lspuce
