    swig -lua -c++ -Iinclude AudioEffectFlanger.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectFlanger.so AudioEffectFlanger_wrap.cxx -lstdc++ -lm -lluajit    
    