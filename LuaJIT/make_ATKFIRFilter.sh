    swig -lua -c++ -Iinclude ATKFIRFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKFIRFilter.so ATKFIRFilter_wrap.cxx -lstdc++ -lm -lluajit    
    