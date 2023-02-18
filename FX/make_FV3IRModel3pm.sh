    swig -lua -c++ -I../include -I../include FV3IRModel3pm.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IRModel3pm.so FV3IRModel3pm_wrap.cxx -lstdc++ -lm -lluajit 
    