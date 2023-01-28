    swig -lua -c++ -Iinclude DelayProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DelayProcessors.so DelayProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    