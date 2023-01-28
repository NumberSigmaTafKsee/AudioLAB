    swig -lua -c++ -Iinclude BlitSaw.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o BlitSaw.so BlitSaw_wrap.cxx -lstdc++ -lm -lluajit    
    