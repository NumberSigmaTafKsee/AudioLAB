    swig -lua -c++ -Iinclude RackFX.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFX.so RackFX_wrap.cxx -lstdc++ -lm -lluajit    
    