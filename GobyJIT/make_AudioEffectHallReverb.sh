    swig -lua -c++ -Iinclude AudioEffectHallReverb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectHallReverb.so AudioEffectHallReverb_wrap.cxx -lstdc++ -lm -lluajit    
    