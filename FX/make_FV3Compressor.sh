    swig -lua -c++ -I../include -I../include FV3Compressor.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3Compressor.so FV3Compressor_wrap.cxx -lstdc++ -lm -lluajit 
    