    swig -lua -c++ -I../include -I../include PeakLimiter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o PeakLimiter.so PeakLimiter_wrap.cxx -lstdc++ -lm -lluajit 
    