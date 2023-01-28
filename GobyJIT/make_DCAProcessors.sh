    swig -lua -c++ -Iinclude DCAProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DCAProcessors.so DCAProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    