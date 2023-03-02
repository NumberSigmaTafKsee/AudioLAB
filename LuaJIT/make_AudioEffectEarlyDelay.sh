    swig -lua -c++ -Iinclude AudioEffectEarlyDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectEarlyDelay.so AudioEffectEarlyDelay_wrap.cxx -lstdc++ -lm -lluajit    
    