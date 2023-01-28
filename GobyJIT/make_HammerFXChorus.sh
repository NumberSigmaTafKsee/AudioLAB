    swig -lua -c++ -Iinclude HammerFXChorus.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXChorus.so HammerFXChorus_wrap.cxx -lstdc++ -lm -lluajit    
    