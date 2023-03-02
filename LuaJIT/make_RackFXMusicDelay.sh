    swig -lua -c++ -Iinclude RackFXMusicDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXMusicDelay.so RackFXMusicDelay_wrap.cxx -lstdc++ -lm -lluajit    
    