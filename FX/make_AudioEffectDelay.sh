    swig -lua -c++ -I../include -I../include AudioEffectDelay.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioEffectDelay.so AudioEffectDelay_wrap.cxx -lstdc++ -lm -lluajit 
    