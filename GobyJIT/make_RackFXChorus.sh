    swig -lua -c++ -Iinclude RackFXChorus.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXChorus.so RackFXChorus_wrap.cxx -lstdc++ -lm -lluajit    
    