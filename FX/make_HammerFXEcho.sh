    swig -lua -c++ -I../include -I../include HammerFXEcho.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o HammerFXEcho.so HammerFXEcho_wrap.cxx -lstdc++ -lm -lluajit 
    