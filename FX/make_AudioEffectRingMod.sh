    swig -lua -c++ -I../include -I../include AudioEffectRingMod.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioEffectRingMod.so AudioEffectRingMod_wrap.cxx -lstdc++ -lm -lluajit 
    