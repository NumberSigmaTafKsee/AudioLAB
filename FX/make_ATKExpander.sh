    swig -lua -c++ -I../include -I../include ATKExpander.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKExpander.so ATKExpander_wrap.cxx -lstdc++ -lm -lluajit 
    