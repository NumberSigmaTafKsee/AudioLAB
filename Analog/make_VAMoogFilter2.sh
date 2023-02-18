    swig -lua -c++ -I../include -I../include/Analog VAMoogFilter2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogFilter2.so VAMoogFilter2_wrap.cxx -lstdc++ -lm -lluajit 
    