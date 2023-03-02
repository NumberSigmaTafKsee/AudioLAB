    swig -lua -c++ -Iinclude Diode.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Diode.so Diode_wrap.cxx -lstdc++ -lm -lluajit    
    