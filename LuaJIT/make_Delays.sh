    swig -lua -c++ -Iinclude Delays.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Delays.so Delays_wrap.cxx -lstdc++ -lm -lluajit    
    