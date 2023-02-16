    swig -lua -c++ -Iinclude VAKorg35LPFFilter.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAKorg35LPFFilter.cpp.so VAKorg35LPFFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    