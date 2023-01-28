    swig -lua -c++ -Iinclude VARKLadderFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VARKLadderFilter.so VARKLadderFilter_wrap.cxx -lstdc++ -lm -lluajit    
    