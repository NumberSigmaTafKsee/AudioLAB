    swig -lua -c++ -Iinclude DelayLine.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DelayLine.so DelayLine_wrap.cxx -lstdc++ -lm -lluajit    
    