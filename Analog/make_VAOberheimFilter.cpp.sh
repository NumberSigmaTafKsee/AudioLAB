    swig -lua -c++ -I../include -I../include/Analog VAOberheimFilter.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAOberheimFilter.cpp.so VAOberheimFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    