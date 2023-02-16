    swig -lua -c++ -Iinclude ATKFixedDelayLine.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKFixedDelayLine.so ATKFixedDelayLine_wrap.cxx -lstdc++ -lm -lluajit    
    