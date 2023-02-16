    swig -lua -c++ -Iinclude Bezier.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Bezier.so Bezier_wrap.cxx -lstdc++ -lm -lluajit    
    