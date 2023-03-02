    swig -lua -c++ -Iinclude Mu45.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Mu45.so Mu45_wrap.cxx -lstdc++ -lm -lluajit    
    