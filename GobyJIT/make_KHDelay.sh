    swig -lua -c++ -Iinclude KHDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o KHDelay.so KHDelay_wrap.cxx -lstdc++ -lm -lluajit    
    