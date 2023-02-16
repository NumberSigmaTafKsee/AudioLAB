    swig -lua -c++ -Iinclude RackFXVibe.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXVibe.so RackFXVibe_wrap.cxx -lstdc++ -lm -lluajit    
    