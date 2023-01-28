swig -lua -c++ spucelpq.i
gcc -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o spucelpq.so spucelpq_wrap.cxx -lstdc++ -lm -lluajit -lspuce
