    swig -lua -c++ -Iinclude Modulation.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Modulation.so Modulation_wrap.cxx -lstdc++ -lm -lluajit    
    