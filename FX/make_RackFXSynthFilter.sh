    swig -lua -c++ -I../include -I../include RackFXSynthFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXSynthFilter.so RackFXSynthFilter_wrap.cxx -lstdc++ -lm -lluajit 
    