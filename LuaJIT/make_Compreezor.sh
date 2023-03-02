    swig -lua -c++ -Iinclude Compreezor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Compreezor.so Compreezor_wrap.cxx -lstdc++ -lm -lluajit    
    