    swig -lua -c++ -Iinclude HammerFXDistortion.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o HammerFXDistortion.so HammerFXDistortion_wrap.cxx -lstdc++ -lm -lluajit    
    