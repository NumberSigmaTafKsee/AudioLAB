    swig -lua -c++ -I../include -I../include HammerFXVibrato.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o HammerFXVibrato.so HammerFXVibrato_wrap.cxx -lstdc++ -lm -lluajit 
    