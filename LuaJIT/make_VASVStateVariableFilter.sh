    swig -lua -c++ -Iinclude VASVStateVariableFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VASVStateVariableFilter.so VASVStateVariableFilter_wrap.cxx -lstdc++ -lm -lluajit    
    