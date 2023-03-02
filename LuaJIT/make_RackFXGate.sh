    swig -lua -c++ -Iinclude RackFXGate.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXGate.so RackFXGate_wrap.cxx -lstdc++ -lm -lluajit    
    