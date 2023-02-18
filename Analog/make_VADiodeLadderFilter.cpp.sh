    swig -lua -c++ -I../include -I../include/Analog VADiodeLadderFilter.cpp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VADiodeLadderFilter.cpp.so VADiodeLadderFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit 
    