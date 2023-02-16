    swig -lua -c++ -Iinclude VAMoogLadders.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogLadders.so VAMoogLadders_wrap.cxx -lstdc++ -lm -lluajit    
    