    swig -lua -c++ -I../include -I../include/Analog VAMoogHalfLadderFilter.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogHalfLadderFilter.cpp.so VAMoogHalfLadderFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    