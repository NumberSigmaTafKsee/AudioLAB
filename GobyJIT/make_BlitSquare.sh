    swig -lua -c++ -Iinclude BlitSquare.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o BlitSquare.so BlitSquare_wrap.cxx -lstdc++ -lm -lluajit    
    