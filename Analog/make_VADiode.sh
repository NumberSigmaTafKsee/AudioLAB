    swig -lua -c++ -I../include -I../include/Analog VADiode.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VADiode.so VADiode_wrap.cxx -lstdc++ -lm -lluajit 
    