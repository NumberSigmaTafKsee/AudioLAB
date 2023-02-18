    swig -lua -c++ -I../include -I../include FV3NRev.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3NRev.so FV3NRev_wrap.cxx -lstdc++ -lm -lluajit 
    