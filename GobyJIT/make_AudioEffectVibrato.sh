    swig -lua -c++ -Iinclude AudioEffectVibrato.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectVibrato.so AudioEffectVibrato_wrap.cxx -lstdc++ -lm -lluajit    
    