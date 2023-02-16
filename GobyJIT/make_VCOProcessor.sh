    swig -lua -c++ -Iinclude VCOProcessor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VCOProcessor.so VCOProcessor_wrap.cxx -lstdc++ -lm -lluajit    
    