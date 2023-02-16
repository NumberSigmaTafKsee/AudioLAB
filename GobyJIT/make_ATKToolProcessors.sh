    swig -lua -c++ -Iinclude ATKToolProcessors.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKToolProcessors.so ATKToolProcessors_wrap.cxx -lstdc++ -lm -lluajit    
    