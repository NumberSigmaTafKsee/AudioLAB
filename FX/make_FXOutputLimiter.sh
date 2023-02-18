    swig -lua -c++ -I../include -I../include FXOutputLimiter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FXOutputLimiter.so FXOutputLimiter_wrap.cxx -lstdc++ -lm -lluajit 
    