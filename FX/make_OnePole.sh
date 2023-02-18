    swig -lua -c++ -I../include -I../include OnePole.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o OnePole.so OnePole_wrap.cxx -lstdc++ -lm -lluajit 
    