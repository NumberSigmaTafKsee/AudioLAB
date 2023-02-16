    swig -lua -c++ -Iinclude ATKGainLimiter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKGainLimiter.so ATKGainLimiter_wrap.cxx -lstdc++ -lm -lluajit    
    