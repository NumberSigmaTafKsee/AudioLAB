    swig -lua -c++ -Iinclude VAMoogLadderFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogLadderFilter.so VAMoogLadderFilter_wrap.cxx -lstdc++ -lm -lluajit    
    