    swig -lua -c++ -I../include -I../include RackFXCompressor.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXCompressor.so RackFXCompressor_wrap.cxx -lstdc++ -lm -lluajit 
    