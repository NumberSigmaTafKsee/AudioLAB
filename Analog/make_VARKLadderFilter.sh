    swig -lua -c++ -I../include -I../include/Analog VARKLadderFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VARKLadderFilter.so VARKLadderFilter_wrap.cxx -lstdc++ -lm -lluajit 
    