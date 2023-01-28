    swig -lua -c++ -Iinclude AudioTK.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioTK.so AudioTK_wrap.cxx -lstdc++ -lm -lluajit    
    