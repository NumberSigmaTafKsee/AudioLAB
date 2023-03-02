    swig -lua -c++ -Iinclude VAVoltageControlledFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAVoltageControlledFilter.so VAVoltageControlledFilter_wrap.cxx -lstdc++ -lm -lluajit    
    