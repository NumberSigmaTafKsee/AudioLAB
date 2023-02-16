    swig -lua -c++ -Iinclude FilterBank.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FilterBank.cpp.so FilterBank.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    