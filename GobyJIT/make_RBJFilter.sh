    swig -lua -c++ -Iinclude RBJFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RBJFilter.so RBJFilter_wrap.cxx -lstdc++ -lm -lluajit    
    