    swig -lua -c++ -Iinclude AudioEffectCompressor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioEffectCompressor.so AudioEffectCompressor_wrap.cxx -lstdc++ -lm -lluajit    
    