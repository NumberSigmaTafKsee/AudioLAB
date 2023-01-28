    swig -lua -c++ -Iinclude DelaySmooth.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DelaySmooth.so DelaySmooth_wrap.cxx -lstdc++ -lm -lluajit    
    