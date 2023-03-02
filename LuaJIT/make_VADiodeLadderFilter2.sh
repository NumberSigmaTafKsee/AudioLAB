    swig -lua -c++ -Iinclude VADiodeLadderFilter2.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VADiodeLadderFilter2.so VADiodeLadderFilter2_wrap.cxx -lstdc++ -lm -lluajit    
    