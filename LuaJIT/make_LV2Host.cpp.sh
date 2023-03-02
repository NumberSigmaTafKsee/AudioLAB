    swig -lua -c++ -Iinclude LV2Host.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o LV2Host.cpp.so LV2Host.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    