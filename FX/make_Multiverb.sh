    swig -lua -c++ -I../include -I../include Multiverb.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Multiverb.so Multiverb_wrap.cxx -lstdc++ -lm -lluajit 
    