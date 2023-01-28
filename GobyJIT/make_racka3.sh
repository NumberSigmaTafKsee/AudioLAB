    swig -lua -c++ -Iinclude racka3.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o racka3.so racka3_wrap.cxx -lstdc++ -lm -lluajit    
    