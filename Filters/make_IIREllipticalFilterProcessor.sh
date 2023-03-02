    swig -lua -c++ -I../include -I../include IIREllipticalFilterProcessor.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o IIREllipticalFilterProcessor.so IIREllipticalFilterProcessor_wrap.cxx -lstdc++ -lm -lluajit 
    