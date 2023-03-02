    swig -lua -c++ -Iinclude HammerFXVibrato.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXVibrato.so HammerFXVibrato_wrap.cxx -lstdc++ -lm -lluajit    
    