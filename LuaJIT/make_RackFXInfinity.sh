    swig -lua -c++ -Iinclude RackFXInfinity.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXInfinity.so RackFXInfinity_wrap.cxx -lstdc++ -lm -lluajit    
    