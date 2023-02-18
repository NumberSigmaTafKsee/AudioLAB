    swig -lua -c++ -I../include -I../include FV3AllPassDecay.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3AllPassDecay.so FV3AllPassDecay_wrap.cxx -lstdc++ -lm -lluajit 
    