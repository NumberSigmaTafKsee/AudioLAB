    swig -lua -c++ -Iinclude Limiter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Limiter.so Limiter_wrap.cxx -lstdc++ -lm -lluajit    
    