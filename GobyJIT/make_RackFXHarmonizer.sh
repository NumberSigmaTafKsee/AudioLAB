    swig -lua -c++ -Iinclude RackFXHarmonizer.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXHarmonizer.so RackFXHarmonizer_wrap.cxx -lstdc++ -lm -lluajit    
    