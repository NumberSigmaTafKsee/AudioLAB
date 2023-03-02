    swig -lua -c++ -Iinclude ATKGainFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKGainFilter.so ATKGainFilter_wrap.cxx -lstdc++ -lm -lluajit    
    