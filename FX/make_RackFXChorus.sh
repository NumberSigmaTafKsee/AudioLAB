    swig -lua -c++ -I../include -I../include RackFXChorus.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXChorus.so RackFXChorus_wrap.cxx -lstdc++ -lm -lluajit 
    