    swig -lua -c++ -Iinclude VAMoogHalfLadderFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogHalfLadderFilter.so VAMoogHalfLadderFilter_wrap.cxx -lstdc++ -lm -lluajit    
    