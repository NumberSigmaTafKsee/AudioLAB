    swig -lua -c++ -I../include -I../include/Analog VAMoogFilter3.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogFilter3.so VAMoogFilter3_wrap.cxx -lstdc++ -lm -lluajit 
    