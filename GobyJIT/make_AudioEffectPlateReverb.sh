    swig -lua -c++ -Iinclude AudioEffectPlateReverb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectPlateReverb.so AudioEffectPlateReverb_wrap.cxx -lstdc++ -lm -lluajit    
    