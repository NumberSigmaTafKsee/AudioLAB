swig -lua -c++ iirspucecheby2bp.i
gcc -O2 -fPIC -march=native -mavx2 -shared -o spucecheby2bp.so iirspucecheby2bp_wrap.cxx -lstdc++ -lm -lluajit -lspuce
