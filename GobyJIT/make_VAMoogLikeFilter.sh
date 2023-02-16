    swig -lua -c++ -Iinclude VAMoogLikeFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogLikeFilter.so VAMoogLikeFilter_wrap.cxx -lstdc++ -lm -lluajit    
    