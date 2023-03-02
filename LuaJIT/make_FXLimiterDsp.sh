    swig -lua -c++ -Iinclude FXLimiterDsp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXLimiterDsp.so FXLimiterDsp_wrap.cxx -lstdc++ -lm -lluajit    
    