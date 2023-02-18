    swig -lua -c++ -I../include -I../include AudioEffectSVF.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioEffectSVF.so AudioEffectSVF_wrap.cxx -lstdc++ -lm -lluajit 
    