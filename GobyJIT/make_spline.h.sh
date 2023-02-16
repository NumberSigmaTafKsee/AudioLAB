    swig -lua -c++ -Iinclude spline.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o spline.h.so spline.h_wrap.cxx -lstdc++ -lm -lluajit    
    