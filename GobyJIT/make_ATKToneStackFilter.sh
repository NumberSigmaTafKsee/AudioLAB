    swig -lua -c++ -Iinclude ATKToneStackFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKToneStackFilter.so ATKToneStackFilter_wrap.cxx -lstdc++ -lm -lluajit    
    