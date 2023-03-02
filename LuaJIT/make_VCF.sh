    swig -lua -c++ -Iinclude VCF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VCF.so VCF_wrap.cxx -lstdc++ -lm -lluajit    
    