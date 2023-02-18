    swig -lua -c++ -I../include -I../include/Analog VAMoogLadderFilter.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogLadderFilter.cpp.so VAMoogLadderFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    