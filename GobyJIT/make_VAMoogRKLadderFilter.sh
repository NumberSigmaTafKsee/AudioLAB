    swig -lua -c++ -Iinclude VAMoogRKLadderFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogRKLadderFilter.so VAMoogRKLadderFilter_wrap.cxx -lstdc++ -lm -lluajit    
    