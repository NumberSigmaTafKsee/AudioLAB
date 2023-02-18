    swig -lua -c++ -I../include -I../include FV3PlainDelay.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3PlainDelay.so FV3PlainDelay_wrap.cxx -lstdc++ -lm -lluajit 
    