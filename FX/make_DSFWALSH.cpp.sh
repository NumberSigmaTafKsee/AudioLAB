    swig -lua -c++ -I../include -I../include DSFWALSH.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o DSFWALSH.cpp.so DSFWALSH.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    