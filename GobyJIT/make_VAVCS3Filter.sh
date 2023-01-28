    swig -lua -c++ -Iinclude VAVCS3Filter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAVCS3Filter.so VAVCS3Filter_wrap.cxx -lstdc++ -lm -lluajit    
    