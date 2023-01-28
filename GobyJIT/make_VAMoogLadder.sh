    swig -lua -c++ -Iinclude VAMoogLadder.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogLadder.so VAMoogLadder_wrap.cxx -lstdc++ -lm -lluajit    
    