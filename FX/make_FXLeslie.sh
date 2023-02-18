    swig -lua -c++ -I../include -I../include FXLeslie.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FXLeslie.so FXLeslie_wrap.cxx -lstdc++ -lm -lluajit 
    