    swig -lua -c++ -Iinclude VAOberheimFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAOberheimFilter.so VAOberheimFilter_wrap.cxx -lstdc++ -lm -lluajit    
    