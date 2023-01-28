    swig -lua -c++ -Iinclude RackFXMBDist.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXMBDist.so RackFXMBDist_wrap.cxx -lstdc++ -lm -lluajit    
    