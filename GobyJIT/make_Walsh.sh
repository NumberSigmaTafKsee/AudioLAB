    swig -lua -c++ -Iinclude Walsh.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Walsh.so Walsh_wrap.cxx -lstdc++ -lm -lluajit    
    