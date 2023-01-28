    swig -lua -c++ -Iinclude Schroeder.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Schroeder.so Schroeder_wrap.cxx -lstdc++ -lm -lluajit    
    