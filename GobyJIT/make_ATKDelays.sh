    swig -lua -c++ -Iinclude ATKDelays.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKDelays.so ATKDelays_wrap.cxx -lstdc++ -lm -lluajit    
    