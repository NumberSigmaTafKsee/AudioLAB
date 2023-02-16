    swig -lua -c++ -Iinclude ATKGainSwell.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKGainSwell.so ATKGainSwell_wrap.cxx -lstdc++ -lm -lluajit    
    