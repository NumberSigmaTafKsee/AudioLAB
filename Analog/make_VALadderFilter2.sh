    swig -lua -c++ -I../include -I../include/Analog VALadderFilter2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VALadderFilter2.so VALadderFilter2_wrap.cxx -lstdc++ -lm -lluajit 
    