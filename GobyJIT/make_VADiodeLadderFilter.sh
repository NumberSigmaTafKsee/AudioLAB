    swig -lua -c++ -Iinclude VADiodeLadderFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VADiodeLadderFilter.so VADiodeLadderFilter_wrap.cxx -lstdc++ -lm -lluajit    
    