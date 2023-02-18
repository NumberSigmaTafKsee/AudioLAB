    swig -lua -c++ -I../include -I../include FV3AllPassModulatedDecay.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3AllPassModulatedDecay.so FV3AllPassModulatedDecay_wrap.cxx -lstdc++ -lm -lluajit 
    