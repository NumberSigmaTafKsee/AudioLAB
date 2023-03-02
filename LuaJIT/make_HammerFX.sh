    swig -lua -c++ -Iinclude HammerFX.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFX.so HammerFX_wrap.cxx -lstdc++ -lm -lluajit    
    