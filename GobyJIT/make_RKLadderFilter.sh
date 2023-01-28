    swig -lua -c++ -Iinclude RKLadderFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RKLadderFilter.so RKLadderFilter_wrap.cxx -lstdc++ -lm -lluajit    
    