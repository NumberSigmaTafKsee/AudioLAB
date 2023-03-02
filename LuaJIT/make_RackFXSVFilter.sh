    swig -lua -c++ -Iinclude RackFXSVFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXSVFilter.so RackFXSVFilter_wrap.cxx -lstdc++ -lm -lluajit    
    