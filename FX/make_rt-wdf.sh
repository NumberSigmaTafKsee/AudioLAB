    swig -lua -c++ -I../include -I../include rt-wdf.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o rt-wdf.so rt-wdf_wrap.cxx -lstdc++ -lm -lluajit 
    