    swig -lua -c++ -I../include -I../include AudioEffectsSuite.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioEffectsSuite.so AudioEffectsSuite_wrap.cxx -lstdc++ -lm -lluajit 
    