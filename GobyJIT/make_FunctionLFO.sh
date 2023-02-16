    swig -lua -c++ -Iinclude FunctionLFO.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FunctionLFO.so FunctionLFO_wrap.cxx -lstdc++ -lm -lluajit    
    