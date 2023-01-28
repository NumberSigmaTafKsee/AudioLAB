    swig -lua -c++ -Iinclude VAMoogHalfLadderFilter.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogHalfLadderFilter.cpp.so VAMoogHalfLadderFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    