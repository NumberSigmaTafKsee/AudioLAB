    swig -lua -c++ -I../include -I../include FV3IRModel2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IRModel2.so FV3IRModel2_wrap.cxx -lstdc++ -lm -lluajit 
    