    swig -lua -c++ -Iinclude RackFXPan.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXPan.so RackFXPan_wrap.cxx -lstdc++ -lm -lluajit    
    