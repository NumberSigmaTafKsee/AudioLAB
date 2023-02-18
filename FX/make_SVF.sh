    swig -lua -c++ -I../include -I../include SVF.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o SVF.so SVF_wrap.cxx -lstdc++ -lm -lluajit 
    