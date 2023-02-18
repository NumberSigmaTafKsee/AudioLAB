    swig -lua -c++ -I../include -I../include RackFXLooper.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXLooper.so RackFXLooper_wrap.cxx -lstdc++ -lm -lluajit 
    