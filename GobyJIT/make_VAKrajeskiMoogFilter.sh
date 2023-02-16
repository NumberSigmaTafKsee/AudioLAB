    swig -lua -c++ -Iinclude VAKrajeskiMoogFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAKrajeskiMoogFilter.so VAKrajeskiMoogFilter_wrap.cxx -lstdc++ -lm -lluajit    
    