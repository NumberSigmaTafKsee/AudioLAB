    swig -lua -c++ -Iinclude VAMoogLadderFilter.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogLadderFilter.cpp.so VAMoogLadderFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    