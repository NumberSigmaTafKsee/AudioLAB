    swig -lua -c++ -Iinclude FXFlanger.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXFlanger.so FXFlanger_wrap.cxx -lstdc++ -lm -lluajit    
    