    swig -lua -c++ -Iinclude ATKAttackReleaseHysterisis.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKAttackReleaseHysterisis.so ATKAttackReleaseHysterisis_wrap.cxx -lstdc++ -lm -lluajit    
    