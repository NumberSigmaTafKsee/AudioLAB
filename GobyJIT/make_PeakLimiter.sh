    swig -lua -c++ -Iinclude PeakLimiter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o PeakLimiter.so PeakLimiter_wrap.cxx -lstdc++ -lm -lluajit    
    