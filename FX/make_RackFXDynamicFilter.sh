    swig -lua -c++ -I../include -I../include RackFXDynamicFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXDynamicFilter.so RackFXDynamicFilter_wrap.cxx -lstdc++ -lm -lluajit 
    