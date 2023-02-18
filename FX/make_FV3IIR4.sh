    swig -lua -c++ -I../include -I../include FV3IIR4.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IIR4.so FV3IIR4_wrap.cxx -lstdc++ -lm -lluajit 
    