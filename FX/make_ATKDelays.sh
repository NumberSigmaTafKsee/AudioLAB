    swig -lua -c++ -I../include -I../include ATKDelays.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKDelays.so ATKDelays_wrap.cxx -lstdc++ -lm -lluajit 
    