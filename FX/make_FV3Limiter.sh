    swig -lua -c++ -I../include -I../include FV3Limiter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3Limiter.so FV3Limiter_wrap.cxx -lstdc++ -lm -lluajit 
    