    swig -lua -c++ -I../include -I../include ATK.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATK.so ATK_wrap.cxx -lstdc++ -lm -lluajit 
    