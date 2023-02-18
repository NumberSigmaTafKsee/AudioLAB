    swig -lua -c++ -I../include -I../include Fir1.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Fir1.cpp.so Fir1.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    