    swig -lua -c++ -Iinclude RackFXCoilCrafter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXCoilCrafter.so RackFXCoilCrafter_wrap.cxx -lstdc++ -lm -lluajit    
    