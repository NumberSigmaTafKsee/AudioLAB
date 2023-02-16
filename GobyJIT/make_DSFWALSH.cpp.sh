    swig -lua -c++ -Iinclude DSFWALSH.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DSFWALSH.cpp.so DSFWALSH.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    