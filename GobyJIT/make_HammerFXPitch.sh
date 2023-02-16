    swig -lua -c++ -Iinclude HammerFXPitch.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXPitch.so HammerFXPitch_wrap.cxx -lstdc++ -lm -lluajit    
    