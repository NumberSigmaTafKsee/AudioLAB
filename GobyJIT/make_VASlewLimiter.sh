    swig -lua -c++ -Iinclude VASlewLimiter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VASlewLimiter.so VASlewLimiter_wrap.cxx -lstdc++ -lm -lluajit    
    