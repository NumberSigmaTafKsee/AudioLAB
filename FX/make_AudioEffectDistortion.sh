    swig -lua -c++ -I../include -I../include AudioEffectDistortion.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AudioEffectDistortion.so AudioEffectDistortion_wrap.cxx -lstdc++ -lm -lluajit 
    