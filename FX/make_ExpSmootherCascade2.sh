    swig -lua -c++ -I../include -I../include ExpSmootherCascade2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ExpSmootherCascade2.so ExpSmootherCascade2_wrap.cxx -lstdc++ -lm -lluajit 
    