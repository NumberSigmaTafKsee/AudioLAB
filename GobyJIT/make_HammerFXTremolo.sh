    swig -lua -c++ -Iinclude HammerFXTremolo.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXTremolo.so HammerFXTremolo_wrap.cxx -lstdc++ -lm -lluajit    
    