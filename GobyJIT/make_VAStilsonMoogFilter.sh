    swig -lua -c++ -Iinclude VAStilsonMoogFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAStilsonMoogFilter.so VAStilsonMoogFilter_wrap.cxx -lstdc++ -lm -lluajit    
    