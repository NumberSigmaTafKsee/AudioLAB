    swig -lua -c++ -Iinclude ATKIIRFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKIIRFilter.so ATKIIRFilter_wrap.cxx -lstdc++ -lm -lluajit    
    