    swig -lua -c++ -I../include -I../include Comb2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Comb2.so Comb2_wrap.cxx -lstdc++ -lm -lluajit 
    