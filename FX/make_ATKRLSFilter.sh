    swig -lua -c++ -I../include -I../include ATKRLSFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKRLSFilter.so ATKRLSFilter_wrap.cxx -lstdc++ -lm -lluajit 
    