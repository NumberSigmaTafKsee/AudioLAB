    swig -lua -c++ -Iinclude VAMoogCatFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogCatFilter.so VAMoogCatFilter_wrap.cxx -lstdc++ -lm -lluajit    
    