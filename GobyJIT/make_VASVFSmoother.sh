    swig -lua -c++ -Iinclude VASVFSmoother.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VASVFSmoother.so VASVFSmoother_wrap.cxx -lstdc++ -lm -lluajit    
    