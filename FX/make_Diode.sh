    swig -lua -c++ -I../include -I../include Diode.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Diode.so Diode_wrap.cxx -lstdc++ -lm -lluajit 
    