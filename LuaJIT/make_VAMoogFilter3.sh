    swig -lua -c++ -Iinclude VAMoogFilter3.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogFilter3.so VAMoogFilter3_wrap.cxx -lstdc++ -lm -lluajit    
    