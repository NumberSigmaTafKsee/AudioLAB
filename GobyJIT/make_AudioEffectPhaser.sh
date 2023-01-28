    swig -lua -c++ -Iinclude AudioEffectPhaser.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectPhaser.so AudioEffectPhaser_wrap.cxx -lstdc++ -lm -lluajit    
    