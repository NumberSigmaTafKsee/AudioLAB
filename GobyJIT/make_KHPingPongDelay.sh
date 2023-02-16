    swig -lua -c++ -Iinclude KHPingPongDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o KHPingPongDelay.so KHPingPongDelay_wrap.cxx -lstdc++ -lm -lluajit    
    