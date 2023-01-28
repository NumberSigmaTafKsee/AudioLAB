    swig -lua -c++ -Iinclude ExpSmootherCascade.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ExpSmootherCascade.so ExpSmootherCascade_wrap.cxx -lstdc++ -lm -lluajit    
    