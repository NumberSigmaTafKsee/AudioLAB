    swig -lua -c++ -Iinclude Fir1.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Fir1.cpp.so Fir1.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    