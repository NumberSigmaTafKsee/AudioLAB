    swig -lua -c++ -Iinclude AudioEffectChorus.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectChorus.so AudioEffectChorus_wrap.cxx -lstdc++ -lm -lluajit    
    