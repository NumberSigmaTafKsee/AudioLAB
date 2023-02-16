    swig -lua -c++ -Iinclude VAMoogNonLinearFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogNonLinearFilter.so VAMoogNonLinearFilter_wrap.cxx -lstdc++ -lm -lluajit    
    