    swig -lua -c++ -I../include -I../include WDFCompressor.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o WDFCompressor.cpp.so WDFCompressor.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    