    swig -lua -c++ -Iinclude OnePole.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o OnePole.so OnePole_wrap.cxx -lstdc++ -lm -lluajit    
    