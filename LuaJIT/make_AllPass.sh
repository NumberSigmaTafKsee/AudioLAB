    swig -lua -c++ -Iinclude AllPass.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AllPass.so AllPass_wrap.cxx -lstdc++ -lm -lluajit    
    