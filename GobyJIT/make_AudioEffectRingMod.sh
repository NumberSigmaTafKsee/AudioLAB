    swig -lua -c++ -Iinclude AudioEffectRingMod.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectRingMod.so AudioEffectRingMod_wrap.cxx -lstdc++ -lm -lluajit    
    