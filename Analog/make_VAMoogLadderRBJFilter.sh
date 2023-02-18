    swig -lua -c++ -I../include -I../include/Analog VAMoogLadderRBJFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogLadderRBJFilter.so VAMoogLadderRBJFilter_wrap.cxx -lstdc++ -lm -lluajit 
    