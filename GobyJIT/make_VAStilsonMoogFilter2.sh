    swig -lua -c++ -Iinclude VAStilsonMoogFilter2.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAStilsonMoogFilter2.so VAStilsonMoogFilter2_wrap.cxx -lstdc++ -lm -lluajit    
    