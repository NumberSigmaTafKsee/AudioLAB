    swig -lua -c++ -I../include -I../include AudioEffectStereoDelay.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioEffectStereoDelay.so AudioEffectStereoDelay_wrap.cxx -lstdc++ -lm -lluajit 
    