    swig -lua -c++ -Iinclude AudioEffectStereoDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectStereoDelay.so AudioEffectStereoDelay_wrap.cxx -lstdc++ -lm -lluajit    
    