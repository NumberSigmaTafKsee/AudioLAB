    swig -lua -c++ -Iinclude VAMoogFilter1.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogFilter1.so VAMoogFilter1_wrap.cxx -lstdc++ -lm -lluajit    
    