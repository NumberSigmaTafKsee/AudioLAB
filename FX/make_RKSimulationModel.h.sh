    swig -lua -c++ -I../include -I../include RKSimulationModel.h.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RKSimulationModel.h.so RKSimulationModel.h_wrap.cxx -lstdc++ -lm -lluajit 
    