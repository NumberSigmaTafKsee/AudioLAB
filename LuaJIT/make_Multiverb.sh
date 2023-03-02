    swig -lua -c++ -Iinclude Multiverb.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Multiverb.so Multiverb_wrap.cxx -lstdc++ -lm -lluajit    
    