    swig -lua -c++ -I../include -I../include DiodeClipper.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o DiodeClipper.so DiodeClipper_wrap.cxx -lstdc++ -lm -lluajit 
    