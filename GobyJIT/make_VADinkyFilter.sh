    swig -lua -c++ -Iinclude VADinkyFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VADinkyFilter.so VADinkyFilter_wrap.cxx -lstdc++ -lm -lluajit    
    