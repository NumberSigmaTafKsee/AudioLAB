    swig -lua -c++ -I../include -I../include Mu45.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Mu45.so Mu45_wrap.cxx -lstdc++ -lm -lluajit 
    