    swig -lua -c++ -Iinclude RackFXWaveshaper.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXWaveshaper.so RackFXWaveshaper_wrap.cxx -lstdc++ -lm -lluajit    
    