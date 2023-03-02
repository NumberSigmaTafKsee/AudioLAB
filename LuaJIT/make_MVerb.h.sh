    swig -lua -c++ -Iinclude MVerb.h.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o MVerb.h.so MVerb.h_wrap.cxx -lstdc++ -lm -lluajit    
    