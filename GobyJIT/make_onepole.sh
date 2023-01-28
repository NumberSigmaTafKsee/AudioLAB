    swig -lua -c++ -Iinclude onepole.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o onepole.so onepole_wrap.cxx -lstdc++ -lm -lluajit    
    