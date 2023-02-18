    swig -lua -c++ -I../include -I../include LFO9000.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o LFO9000.so LFO9000_wrap.cxx -lstdc++ -lm -lluajit 
    