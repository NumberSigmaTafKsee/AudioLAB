    swig -lua -c++ -Iinclude AudioEffectAutoWah.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectAutoWah.so AudioEffectAutoWah_wrap.cxx -lstdc++ -lm -lluajit    
    