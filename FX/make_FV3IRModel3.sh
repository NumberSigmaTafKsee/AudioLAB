    swig -lua -c++ -I../include -I../include FV3IRModel3.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IRModel3.so FV3IRModel3_wrap.cxx -lstdc++ -lm -lluajit 
    