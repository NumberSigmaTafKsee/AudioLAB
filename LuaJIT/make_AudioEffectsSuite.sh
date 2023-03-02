    swig -lua -c++ -Iinclude AudioEffectsSuite.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectsSuite.so AudioEffectsSuite_wrap.cxx -lstdc++ -lm -lluajit    
    