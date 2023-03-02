    swig -lua -c++ -Iinclude RackFXDynamicFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXDynamicFilter.so RackFXDynamicFilter_wrap.cxx -lstdc++ -lm -lluajit    
    