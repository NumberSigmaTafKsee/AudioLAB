    swig -lua -c++ -Iinclude VAOberheimFilter.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAOberheimFilter.cpp.so VAOberheimFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    