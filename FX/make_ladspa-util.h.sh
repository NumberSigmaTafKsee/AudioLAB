    swig -lua -c++ -I../include -I../include ladspa-util.h.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ladspa-util.h.so ladspa-util.h_wrap.cxx -lstdc++ -lm -lluajit 
    