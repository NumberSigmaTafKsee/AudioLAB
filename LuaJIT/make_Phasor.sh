    swig -lua -c++ -Iinclude Phasor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Phasor.so Phasor_wrap.cxx -lstdc++ -lm -lluajit    
    