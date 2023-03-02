    swig -lua -c++ -Iinclude HammerFXRotary.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXRotary.so HammerFXRotary_wrap.cxx -lstdc++ -lm -lluajit    
    