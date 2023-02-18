    swig -lua -c++ -I../include -I../include FXDelays.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FXDelays.so FXDelays_wrap.cxx -lstdc++ -lm -lluajit 
    