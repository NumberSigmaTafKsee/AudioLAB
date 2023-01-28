swig -lua -c++ iirspucebutterlp.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o spucebutterlp.so iirspucebutterlp_wrap.cxx -lstdc++ -lm -lluajit -lspuce
