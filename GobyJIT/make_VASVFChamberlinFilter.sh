    swig -lua -c++ -Iinclude VASVFChamberlinFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VASVFChamberlinFilter.so VASVFChamberlinFilter_wrap.cxx -lstdc++ -lm -lluajit    
    