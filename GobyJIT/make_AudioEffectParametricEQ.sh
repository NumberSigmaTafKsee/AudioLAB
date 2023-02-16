    swig -lua -c++ -Iinclude AudioEffectParametricEQ.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectParametricEQ.so AudioEffectParametricEQ_wrap.cxx -lstdc++ -lm -lluajit    
    