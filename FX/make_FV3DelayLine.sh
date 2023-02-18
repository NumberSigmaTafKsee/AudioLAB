    swig -lua -c++ -I../include -I../include FV3DelayLine.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3DelayLine.so FV3DelayLine_wrap.cxx -lstdc++ -lm -lluajit 
    