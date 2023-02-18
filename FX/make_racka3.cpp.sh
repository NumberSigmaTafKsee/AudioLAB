    swig -lua -c++ -I../include -I../include racka3.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o racka3.cpp.so racka3.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    