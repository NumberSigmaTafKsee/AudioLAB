    swig -lua -c++ -I../include -I../include/Analog VoltageControlledOscillator.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VoltageControlledOscillator.so VoltageControlledOscillator_wrap.cxx -lstdc++ -lm -lluajit 
    