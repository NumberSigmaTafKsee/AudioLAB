    swig -lua -c++ -Iinclude Rubberband.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Rubberband.so Rubberband_wrap.cxx -lstdc++ -lm -lluajit    
    