    swig -lua -c++ -I../include -I../include/Analog VAMoogHalfLadderFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogHalfLadderFilter.so VAMoogHalfLadderFilter_wrap.cxx -lstdc++ -lm -lluajit 
    