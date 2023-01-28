    swig -lua -c++ -Iinclude VALadderFilter2.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VALadderFilter2.so VALadderFilter2_wrap.cxx -lstdc++ -lm -lluajit    
    