    swig -lua -c++ -Iinclude DPW.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DPW.so DPW_wrap.cxx -lstdc++ -lm -lluajit    
    