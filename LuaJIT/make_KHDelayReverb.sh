    swig -lua -c++ -Iinclude KHDelayReverb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o KHDelayReverb.so KHDelayReverb_wrap.cxx -lstdc++ -lm -lluajit    
    