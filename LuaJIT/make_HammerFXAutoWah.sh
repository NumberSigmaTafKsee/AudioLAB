    swig -lua -c++ -Iinclude HammerFXAutoWah.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXAutoWah.so HammerFXAutoWah_wrap.cxx -lstdc++ -lm -lluajit    
    