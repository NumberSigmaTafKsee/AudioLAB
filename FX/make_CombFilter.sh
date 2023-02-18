    swig -lua -c++ -I../include -I../include CombFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o CombFilter.so CombFilter_wrap.cxx -lstdc++ -lm -lluajit 
    