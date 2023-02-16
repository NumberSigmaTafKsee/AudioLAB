    swig -lua -c++ -Iinclude AudioEffectAmplifier.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectAmplifier.so AudioEffectAmplifier_wrap.cxx -lstdc++ -lm -lluajit    
    