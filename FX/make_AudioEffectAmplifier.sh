    swig -lua -c++ -I../include -I../include AudioEffectAmplifier.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioEffectAmplifier.so AudioEffectAmplifier_wrap.cxx -lstdc++ -lm -lluajit 
    