    swig -lua -c++ -Iinclude AudioEffectDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectDelay.so AudioEffectDelay_wrap.cxx -lstdc++ -lm -lluajit    
    