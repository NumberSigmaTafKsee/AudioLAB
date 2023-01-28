    swig -lua -c++ -Iinclude HammerFXEcho.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXEcho.so HammerFXEcho_wrap.cxx -lstdc++ -lm -lluajit    
    