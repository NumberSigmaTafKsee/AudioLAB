    swig -lua -c++ -Iinclude RCFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RCFilter.so RCFilter_wrap.cxx -lstdc++ -lm -lluajit    
    