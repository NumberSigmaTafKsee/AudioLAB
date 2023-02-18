    swig -lua -c++ -I../include -I../include/Analog VAWDFCompressor.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAWDFCompressor.so VAWDFCompressor_wrap.cxx -lstdc++ -lm -lluajit 
    