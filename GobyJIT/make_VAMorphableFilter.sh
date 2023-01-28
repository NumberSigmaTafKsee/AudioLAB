    swig -lua -c++ -Iinclude VAMorphableFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMorphableFilter.so VAMorphableFilter_wrap.cxx -lstdc++ -lm -lluajit    
    