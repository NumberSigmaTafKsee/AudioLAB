    swig -lua -c++ -Iinclude AudioEffectPingPong.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectPingPong.so AudioEffectPingPong_wrap.cxx -lstdc++ -lm -lluajit    
    