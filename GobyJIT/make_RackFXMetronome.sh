    swig -lua -c++ -Iinclude RackFXMetronome.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXMetronome.so RackFXMetronome_wrap.cxx -lstdc++ -lm -lluajit    
    