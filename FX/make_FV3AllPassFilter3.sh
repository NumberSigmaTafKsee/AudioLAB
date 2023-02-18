    swig -lua -c++ -I../include -I../include FV3AllPassFilter3.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3AllPassFilter3.so FV3AllPassFilter3_wrap.cxx -lstdc++ -lm -lluajit 
    