    swig -lua -c++ -Iinclude RackFXLooper.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXLooper.so RackFXLooper_wrap.cxx -lstdc++ -lm -lluajit    
    