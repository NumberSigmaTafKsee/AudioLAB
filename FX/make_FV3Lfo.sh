    swig -lua -c++ -I../include -I../include FV3Lfo.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3Lfo.so FV3Lfo_wrap.cxx -lstdc++ -lm -lluajit 
    