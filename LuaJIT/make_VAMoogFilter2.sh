    swig -lua -c++ -Iinclude VAMoogFilter2.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogFilter2.so VAMoogFilter2_wrap.cxx -lstdc++ -lm -lluajit    
    