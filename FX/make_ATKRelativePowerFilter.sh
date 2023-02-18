    swig -lua -c++ -I../include -I../include ATKRelativePowerFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKRelativePowerFilter.so ATKRelativePowerFilter_wrap.cxx -lstdc++ -lm -lluajit 
    