    swig -lua -c++ -Iinclude RingBuffer.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RingBuffer.so RingBuffer_wrap.cxx -lstdc++ -lm -lluajit    
    