    swig -lua -c++ -Iinclude FXFreeVerb3.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXFreeVerb3.so FXFreeVerb3_wrap.cxx -lstdc++ -lm -lluajit    
    