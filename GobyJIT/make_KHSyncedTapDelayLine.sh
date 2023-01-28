    swig -lua -c++ -Iinclude KHSyncedTapDelayLine.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o KHSyncedTapDelayLine.so KHSyncedTapDelayLine_wrap.cxx -lstdc++ -lm -lluajit    
    