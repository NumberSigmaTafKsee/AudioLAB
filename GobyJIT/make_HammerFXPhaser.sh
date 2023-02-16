    swig -lua -c++ -Iinclude HammerFXPhaser.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXPhaser.so HammerFXPhaser_wrap.cxx -lstdc++ -lm -lluajit    
    