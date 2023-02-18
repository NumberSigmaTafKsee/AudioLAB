    swig -lua -c++ -I../include -I../include/Analog VAKorg35LPFFilter.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAKorg35LPFFilter.cpp.so VAKorg35LPFFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    