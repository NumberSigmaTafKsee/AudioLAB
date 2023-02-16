    swig -lua -c++ -Iinclude ladspa-util.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ladspa-util.h.so ladspa-util.h_wrap.cxx -lstdc++ -lm -lluajit    
    