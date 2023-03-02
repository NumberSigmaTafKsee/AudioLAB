    swig -lua -c++ -Iinclude RackFXEqualizer.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXEqualizer.so RackFXEqualizer_wrap.cxx -lstdc++ -lm -lluajit    
    