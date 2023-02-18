    swig -lua -c++ -I../include -I../include RackFXPan.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXPan.so RackFXPan_wrap.cxx -lstdc++ -lm -lluajit 
    