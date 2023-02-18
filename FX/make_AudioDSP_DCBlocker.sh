    swig -lua -c++ -I../include -I../include AudioDSP_DCBlocker.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioDSP_DCBlocker.so AudioDSP_DCBlocker_wrap.cxx -lstdc++ -lm -lluajit 
    