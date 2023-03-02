    swig -lua -c++ -Iinclude ClippingFunctions.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ClippingFunctions.so ClippingFunctions_wrap.cxx -lstdc++ -lm -lluajit    
    