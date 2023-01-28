    swig -lua -c++ -Iinclude TSTone.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o TSTone.so TSTone_wrap.cxx -lstdc++ -lm -lluajit    
    