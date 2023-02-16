    swig -lua -c++ -Iinclude AudioEffectTremolo.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectTremolo.so AudioEffectTremolo_wrap.cxx -lstdc++ -lm -lluajit    
    