    swig -lua -c++ -Iinclude RKSimulationModel.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RKSimulationModel.h.so RKSimulationModel.h_wrap.cxx -lstdc++ -lm -lluajit    
    