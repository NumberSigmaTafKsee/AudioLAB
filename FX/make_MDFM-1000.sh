    swig -lua -c++ -I../include -I../include MDFM-1000.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o MDFM-1000.so MDFM-1000_wrap.cxx -lstdc++ -lm -lluajit 
    