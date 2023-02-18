    swig -lua -c++ -I../include -I../include/Analog VAStilsonMoogFilter2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAStilsonMoogFilter2.so VAStilsonMoogFilter2_wrap.cxx -lstdc++ -lm -lluajit 
    