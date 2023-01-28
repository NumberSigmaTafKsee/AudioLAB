    swig -lua -c++ -Iinclude DelaySyncedTapDelayLine.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DelaySyncedTapDelayLine.so DelaySyncedTapDelayLine_wrap.cxx -lstdc++ -lm -lluajit    
    