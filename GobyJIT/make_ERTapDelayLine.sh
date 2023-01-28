    swig -lua -c++ -Iinclude ERTapDelayLine.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ERTapDelayLine.so ERTapDelayLine_wrap.cxx -lstdc++ -lm -lluajit    
    