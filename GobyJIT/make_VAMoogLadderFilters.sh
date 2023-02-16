    swig -lua -c++ -Iinclude VAMoogLadderFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogLadderFilters.so VAMoogLadderFilters_wrap.cxx -lstdc++ -lm -lluajit    
    