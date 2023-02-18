    swig -lua -c++ -I../include -I../include FV3IRModel2zm.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IRModel2zm.so FV3IRModel2zm_wrap.cxx -lstdc++ -lm -lluajit 
    