    swig -lua -c++ -I../include -I../include LV2Host.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o LV2Host.so LV2Host_wrap.cxx -lstdc++ -lm -lluajit 
    