    swig -lua -c++ -Iinclude SVF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o SVF.so SVF_wrap.cxx -lstdc++ -lm -lluajit    
    