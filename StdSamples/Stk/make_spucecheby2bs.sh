swig -lua -c++ iirspucecheby2bs.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o spucecheby2bs.so iirspucecheby2bs_wrap.cxx -lstdc++ -lm -lluajit -lspuce
