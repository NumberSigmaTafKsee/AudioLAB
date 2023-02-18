    swig -lua -c++ -I../include -I../include SchroederImproved.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o SchroederImproved.so SchroederImproved_wrap.cxx -lstdc++ -lm -lluajit 
    