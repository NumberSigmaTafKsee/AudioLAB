    swig -lua -c++ -Iinclude AudioEffectEarlyReflection.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectEarlyReflection.so AudioEffectEarlyReflection_wrap.cxx -lstdc++ -lm -lluajit    
    