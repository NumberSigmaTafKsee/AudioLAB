    swig -lua -c++ -Iinclude KHCrossDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o KHCrossDelay.so KHCrossDelay_wrap.cxx -lstdc++ -lm -lluajit    
    