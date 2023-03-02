    swig -lua -c++ -Iinclude HammerFXTubeAmp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXTubeAmp.so HammerFXTubeAmp_wrap.cxx -lstdc++ -lm -lluajit    
    