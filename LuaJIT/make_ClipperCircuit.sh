    swig -lua -c++ -Iinclude ClipperCircuit.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ClipperCircuit.so ClipperCircuit_wrap.cxx -lstdc++ -lm -lluajit    
    