    swig -lua -c++ -I../include -I../include RackFXInfinity.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXInfinity.so RackFXInfinity_wrap.cxx -lstdc++ -lm -lluajit 
    