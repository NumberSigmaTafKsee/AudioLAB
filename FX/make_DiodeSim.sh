    swig -lua -c++ -I../include -I../include DiodeSim.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o DiodeSim.so DiodeSim_wrap.cxx -lstdc++ -lm -lluajit 
    