    swig -lua -c++ -Iinclude AudioEffectLFO.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectLFO.so AudioEffectLFO_wrap.cxx -lstdc++ -lm -lluajit    
    