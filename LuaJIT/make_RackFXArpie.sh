    swig -lua -c++ -Iinclude RackFXArpie.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXArpie.so RackFXArpie_wrap.cxx -lstdc++ -lm -lluajit    
    