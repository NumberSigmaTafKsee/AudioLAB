    swig -lua -c++ -Iinclude RackFXAlienWah.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXAlienWah.so RackFXAlienWah_wrap.cxx -lstdc++ -lm -lluajit    
    