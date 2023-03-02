    swig -lua -c++ -Iinclude ATKMultipleUniversalFixedDelayLine.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKMultipleUniversalFixedDelayLine.so ATKMultipleUniversalFixedDelayLine_wrap.cxx -lstdc++ -lm -lluajit    
    