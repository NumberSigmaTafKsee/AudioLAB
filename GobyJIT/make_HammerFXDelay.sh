    swig -lua -c++ -Iinclude HammerFXDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXDelay.so HammerFXDelay_wrap.cxx -lstdc++ -lm -lluajit    
    