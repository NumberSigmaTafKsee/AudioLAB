    swig -lua -c++ -I../include -I../include FXDynamics.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FXDynamics.so FXDynamics_wrap.cxx -lstdc++ -lm -lluajit 
    