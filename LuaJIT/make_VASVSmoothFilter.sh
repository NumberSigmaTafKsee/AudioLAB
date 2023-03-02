    swig -lua -c++ -Iinclude VASVSmoothFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VASVSmoothFilter.so VASVSmoothFilter_wrap.cxx -lstdc++ -lm -lluajit    
    