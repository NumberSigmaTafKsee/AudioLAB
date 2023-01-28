    swig -lua -c++ -Iinclude ATKPowerFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKPowerFilter.so ATKPowerFilter_wrap.cxx -lstdc++ -lm -lluajit    
    