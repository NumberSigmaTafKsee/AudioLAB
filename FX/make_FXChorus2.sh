    swig -lua -c++ -I../include -I../include FXChorus2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FXChorus2.so FXChorus2_wrap.cxx -lstdc++ -lm -lluajit 
    