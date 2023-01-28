    swig -lua -c++ -Iinclude RingBuffers.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RingBuffers.so RingBuffers_wrap.cxx -lstdc++ -lm -lluajit    
    