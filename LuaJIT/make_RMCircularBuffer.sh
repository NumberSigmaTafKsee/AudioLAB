    swig -lua -c++ -Iinclude RMCircularBuffer.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RMCircularBuffer.so RMCircularBuffer_wrap.cxx -lstdc++ -lm -lluajit    
    