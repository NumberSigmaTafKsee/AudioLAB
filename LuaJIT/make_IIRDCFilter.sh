    swig -lua -c++ -Iinclude IIRDCFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o IIRDCFilter.so IIRDCFilter_wrap.cxx -lstdc++ -lm -lluajit    
    