    swig -lua -c++ -Iinclude VAMoogNonLinearFilter2.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAMoogNonLinearFilter2.so VAMoogNonLinearFilter2_wrap.cxx -lstdc++ -lm -lluajit    
    