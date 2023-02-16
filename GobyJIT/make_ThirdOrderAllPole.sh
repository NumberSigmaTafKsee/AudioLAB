    swig -lua -c++ -Iinclude ThirdOrderAllPole.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ThirdOrderAllPole.so ThirdOrderAllPole_wrap.cxx -lstdc++ -lm -lluajit    
    