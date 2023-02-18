    swig -lua -c++ -I../include -I../include FV3CombFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3CombFilter.so FV3CombFilter_wrap.cxx -lstdc++ -lm -lluajit 
    