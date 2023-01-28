    swig -lua -c++ -Iinclude RackFXResample.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXResample.so RackFXResample_wrap.cxx -lstdc++ -lm -lluajit    
    