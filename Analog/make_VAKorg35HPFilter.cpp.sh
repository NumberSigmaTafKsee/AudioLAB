    swig -lua -c++ -I../include -I../include/Analog VAKorg35HPFilter.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAKorg35HPFilter.cpp.so VAKorg35HPFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    