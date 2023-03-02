    swig -lua -c++ -Iinclude VAMoogFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogFilter.so VAMoogFilter_wrap.cxx -lstdc++ -lm -lluajit    
    