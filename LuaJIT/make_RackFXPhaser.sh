    swig -lua -c++ -Iinclude RackFXPhaser.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXPhaser.so RackFXPhaser_wrap.cxx -lstdc++ -lm -lluajit    
    