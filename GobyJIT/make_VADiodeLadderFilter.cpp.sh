    swig -lua -c++ -Iinclude VADiodeLadderFilter.cpp.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VADiodeLadderFilter.cpp.so VADiodeLadderFilter.cpp_wrap.cxx -lstdc++ -lm -lluajit    
    