    swig -lua -c++ -Iinclude rt-wdf.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o rt-wdf.so rt-wdf_wrap.cxx -lstdc++ -lm -lluajit    
    