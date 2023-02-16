    swig -lua -c++ -Iinclude RackFXDFlange.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXDFlange.so RackFXDFlange_wrap.cxx -lstdc++ -lm -lluajit    
    