    swig -lua -c++ -I../include -I../include FXOscillators.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FXOscillators.so FXOscillators_wrap.cxx -lstdc++ -lm -lluajit 
    