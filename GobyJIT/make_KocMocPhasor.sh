    swig -lua -c++ -Iinclude KocMocPhasor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o KocMocPhasor.so KocMocPhasor_wrap.cxx -lstdc++ -lm -lluajit    
    