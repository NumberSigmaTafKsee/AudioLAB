    swig -lua -c++ -Iinclude ..i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ..so ._wrap.cxx -lstdc++ -lm -lluajit    
    