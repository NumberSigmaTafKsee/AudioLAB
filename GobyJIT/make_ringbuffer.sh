    swig -lua -c++ -Iinclude ringbuffer.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ringbuffer.so ringbuffer_wrap.cxx -lstdc++ -lm -lluajit    
    