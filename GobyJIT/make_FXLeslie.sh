    swig -lua -c++ -Iinclude FXLeslie.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXLeslie.so FXLeslie_wrap.cxx -lstdc++ -lm -lluajit    
    