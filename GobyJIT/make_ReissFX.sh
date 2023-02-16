    swig -lua -c++ -Iinclude ReissFX.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ReissFX.so ReissFX_wrap.cxx -lstdc++ -lm -lluajit    
    