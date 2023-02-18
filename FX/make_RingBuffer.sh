    swig -lua -c++ -I../include -I../include RingBuffer.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RingBuffer.so RingBuffer_wrap.cxx -lstdc++ -lm -lluajit 
    