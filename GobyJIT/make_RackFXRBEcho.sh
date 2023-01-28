    swig -lua -c++ -Iinclude RackFXRBEcho.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXRBEcho.so RackFXRBEcho_wrap.cxx -lstdc++ -lm -lluajit    
    