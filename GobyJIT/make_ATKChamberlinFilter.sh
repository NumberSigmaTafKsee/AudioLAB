    swig -lua -c++ -Iinclude ATKChamberlinFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKChamberlinFilter.so ATKChamberlinFilter_wrap.cxx -lstdc++ -lm -lluajit    
    