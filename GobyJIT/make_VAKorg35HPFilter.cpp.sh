    swig -lua -c++ -Iinclude VAKorg35HPFilter.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAKorg35HPFilter.cpp.so VAKorg35HPFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    