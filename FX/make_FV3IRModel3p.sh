    swig -lua -c++ -I../include -I../include FV3IRModel3p.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IRModel3p.so FV3IRModel3p_wrap.cxx -lstdc++ -lm -lluajit 
    