    swig -lua -c++ -Iinclude racka3.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o racka3.cpp.so racka3.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    