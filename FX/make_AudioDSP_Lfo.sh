    swig -lua -c++ -I../include -I../include AudioDSP_Lfo.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioDSP_Lfo.so AudioDSP_Lfo_wrap.cxx -lstdc++ -lm -lluajit 
    