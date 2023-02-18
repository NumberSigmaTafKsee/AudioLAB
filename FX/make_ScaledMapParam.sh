    swig -lua -c++ -I../include -I../include ScaledMapParam.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ScaledMapParam.so ScaledMapParam_wrap.cxx -lstdc++ -lm -lluajit 
    