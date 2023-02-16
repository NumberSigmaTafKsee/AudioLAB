    swig -lua -c++ -Iinclude FXOutputLimiter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXOutputLimiter.so FXOutputLimiter_wrap.cxx -lstdc++ -lm -lluajit    
    