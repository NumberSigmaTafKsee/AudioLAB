    swig -lua -c++ -Iinclude RingBufferProcessor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RingBufferProcessor.so RingBufferProcessor_wrap.cxx -lstdc++ -lm -lluajit    
    