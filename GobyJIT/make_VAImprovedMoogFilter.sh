    swig -lua -c++ -Iinclude VAImprovedMoogFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAImprovedMoogFilter.so VAImprovedMoogFilter_wrap.cxx -lstdc++ -lm -lluajit    
    