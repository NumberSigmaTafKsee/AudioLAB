    swig -lua -c++ -Iinclude ATKAttackRelease.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKAttackRelease.so ATKAttackRelease_wrap.cxx -lstdc++ -lm -lluajit    
    