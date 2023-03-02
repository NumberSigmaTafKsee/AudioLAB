    swig -lua -c++ -Iinclude DelayLines.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DelayLines.so DelayLines_wrap.cxx -lstdc++ -lm -lluajit    
    