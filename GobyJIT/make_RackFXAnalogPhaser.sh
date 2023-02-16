    swig -lua -c++ -Iinclude RackFXAnalogPhaser.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXAnalogPhaser.so RackFXAnalogPhaser_wrap.cxx -lstdc++ -lm -lluajit    
    