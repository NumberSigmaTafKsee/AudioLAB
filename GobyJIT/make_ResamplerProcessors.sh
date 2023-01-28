    swig -lua -c++ -Iinclude ResamplerProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ResamplerProcessors.so ResamplerProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    