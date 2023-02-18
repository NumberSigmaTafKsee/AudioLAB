    swig -lua -c++ -I../include -I../include FV3SLimiter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3SLimiter.so FV3SLimiter_wrap.cxx -lstdc++ -lm -lluajit 
    