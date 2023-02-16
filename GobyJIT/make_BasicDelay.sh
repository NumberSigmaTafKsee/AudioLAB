    swig -lua -c++ -Iinclude BasicDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o BasicDelay.so BasicDelay_wrap.cxx -lstdc++ -lm -lluajit    
    