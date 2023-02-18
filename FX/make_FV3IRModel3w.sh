    swig -lua -c++ -I../include -I../include FV3IRModel3w.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IRModel3w.so FV3IRModel3w_wrap.cxx -lstdc++ -lm -lluajit 
    