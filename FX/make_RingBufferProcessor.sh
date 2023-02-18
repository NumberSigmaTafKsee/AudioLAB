    swig -lua -c++ -I../include -I../include RingBufferProcessor.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RingBufferProcessor.so RingBufferProcessor_wrap.cxx -lstdc++ -lm -lluajit 
    