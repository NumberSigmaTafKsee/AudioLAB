    swig -lua -c++ -Iinclude AudioEffectSVF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectSVF.so AudioEffectSVF_wrap.cxx -lstdc++ -lm -lluajit    
    