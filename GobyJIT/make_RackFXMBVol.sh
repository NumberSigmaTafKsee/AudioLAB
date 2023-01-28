    swig -lua -c++ -Iinclude RackFXMBVol.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXMBVol.so RackFXMBVol_wrap.cxx -lstdc++ -lm -lluajit    
    