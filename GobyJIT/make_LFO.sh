    swig -lua -c++ -Iinclude LFO.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o LFO.so LFO_wrap.cxx -lstdc++ -lm -lluajit    
    