    swig -lua -c++ -Iinclude RackFXSynthFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXSynthFilter.so RackFXSynthFilter_wrap.cxx -lstdc++ -lm -lluajit    
    