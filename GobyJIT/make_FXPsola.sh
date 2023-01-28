    swig -lua -c++ -Iinclude FXPsola.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXPsola.so FXPsola_wrap.cxx -lstdc++ -lm -lluajit    
    