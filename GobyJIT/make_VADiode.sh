    swig -lua -c++ -Iinclude VADiode.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VADiode.so VADiode_wrap.cxx -lstdc++ -lm -lluajit    
    