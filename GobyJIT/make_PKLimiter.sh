    swig -lua -c++ -Iinclude PKLimiter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o PKLimiter.so PKLimiter_wrap.cxx -lstdc++ -lm -lluajit    
    