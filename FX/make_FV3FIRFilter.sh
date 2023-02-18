    swig -lua -c++ -I../include -I../include FV3FIRFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3FIRFilter.so FV3FIRFilter_wrap.cxx -lstdc++ -lm -lluajit 
    