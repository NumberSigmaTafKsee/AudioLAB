    swig -lua -c++ -Iinclude VCFProcessor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VCFProcessor.so VCFProcessor_wrap.cxx -lstdc++ -lm -lluajit    
    