    swig -lua -c++ -I../include -I../include AudioDSP_FirstOrderAllPass.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioDSP_FirstOrderAllPass.so AudioDSP_FirstOrderAllPass_wrap.cxx -lstdc++ -lm -lluajit 
    