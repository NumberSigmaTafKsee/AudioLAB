    swig -lua -c++ -Iinclude Moorer.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Moorer.so Moorer_wrap.cxx -lstdc++ -lm -lluajit    
    