    swig -lua -c++ -I../include -I../include AudioDSP_Delay.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioDSP_Delay.so AudioDSP_Delay_wrap.cxx -lstdc++ -lm -lluajit 
    