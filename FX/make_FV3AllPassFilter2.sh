    swig -lua -c++ -I../include -I../include FV3AllPassFilter2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3AllPassFilter2.so FV3AllPassFilter2_wrap.cxx -lstdc++ -lm -lluajit 
    