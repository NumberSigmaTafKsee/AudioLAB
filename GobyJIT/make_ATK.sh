    swig -lua -c++ -Iinclude ATK.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATK.so ATK_wrap.cxx -lstdc++ -lm -lluajit    
    