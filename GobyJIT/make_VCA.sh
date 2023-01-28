    swig -lua -c++ -Iinclude VCA.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VCA.so VCA_wrap.cxx -lstdc++ -lm -lluajit    
    