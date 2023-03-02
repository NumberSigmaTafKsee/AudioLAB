    swig -lua -c++ -Iinclude LFO9000.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o LFO9000.so LFO9000_wrap.cxx -lstdc++ -lm -lluajit    
    