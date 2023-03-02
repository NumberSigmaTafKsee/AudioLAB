    swig -lua -c++ -Iinclude RBJ.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RBJ.so RBJ_wrap.cxx -lstdc++ -lm -lluajit    
    