    swig -lua -c++ -Iinclude FunctionGenerator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FunctionGenerator.so FunctionGenerator_wrap.cxx -lstdc++ -lm -lluajit    
    