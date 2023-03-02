    swig -lua -c++ -Iinclude SlewLimiter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o SlewLimiter.so SlewLimiter_wrap.cxx -lstdc++ -lm -lluajit    
    