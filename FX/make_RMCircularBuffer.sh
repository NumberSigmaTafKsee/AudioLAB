    swig -lua -c++ -I../include -I../include RMCircularBuffer.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RMCircularBuffer.so RMCircularBuffer_wrap.cxx -lstdc++ -lm -lluajit 
    