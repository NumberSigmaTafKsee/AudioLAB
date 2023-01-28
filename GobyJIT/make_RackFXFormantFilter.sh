    swig -lua -c++ -Iinclude RackFXFormantFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXFormantFilter.so RackFXFormantFilter_wrap.cxx -lstdc++ -lm -lluajit    
    