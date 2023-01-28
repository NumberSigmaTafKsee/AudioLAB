    swig -lua -c++ -Iinclude Comb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Comb.so Comb_wrap.cxx -lstdc++ -lm -lluajit    
    