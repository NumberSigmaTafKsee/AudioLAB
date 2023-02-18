    swig -lua -c++ -I../include -I../include FV3IIR1.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IIR1.so FV3IIR1_wrap.cxx -lstdc++ -lm -lluajit 
    