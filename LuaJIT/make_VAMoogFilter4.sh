    swig -lua -c++ -Iinclude VAMoogFilter4.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogFilter4.so VAMoogFilter4_wrap.cxx -lstdc++ -lm -lluajit    
    