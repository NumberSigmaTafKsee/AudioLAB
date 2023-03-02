    swig -lua -c++ -Iinclude ATKExpander.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKExpander.so ATKExpander_wrap.cxx -lstdc++ -lm -lluajit    
    