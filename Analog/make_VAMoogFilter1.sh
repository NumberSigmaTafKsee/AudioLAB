    swig -lua -c++ -I../include -I../include/Analog VAMoogFilter1.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogFilter1.so VAMoogFilter1_wrap.cxx -lstdc++ -lm -lluajit 
    