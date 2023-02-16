    swig -lua -c++ -Iinclude RackFXRBFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o RackFXRBFilter.so RackFXRBFilter_wrap.cxx -lstdc++ -lm -lluajit    
    