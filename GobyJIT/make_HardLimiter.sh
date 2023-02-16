    swig -lua -c++ -Iinclude HardLimiter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HardLimiter.so HardLimiter_wrap.cxx -lstdc++ -lm -lluajit    
    