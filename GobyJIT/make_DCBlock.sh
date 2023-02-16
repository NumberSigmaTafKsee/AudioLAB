    swig -lua -c++ -Iinclude DCBlock.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DCBlock.so DCBlock_wrap.cxx -lstdc++ -lm -lluajit    
    