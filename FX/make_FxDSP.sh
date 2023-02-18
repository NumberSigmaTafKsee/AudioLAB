    swig -lua -c++ -I../include -I../include FxDSP.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FxDSP.so FxDSP_wrap.cxx -lstdc++ -lm -lluajit 
    