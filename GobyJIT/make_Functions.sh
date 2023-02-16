    swig -lua -c++ -Iinclude Functions.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Functions.so Functions_wrap.cxx -lstdc++ -lm -lluajit    
    