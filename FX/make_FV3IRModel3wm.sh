    swig -lua -c++ -I../include -I../include FV3IRModel3wm.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IRModel3wm.so FV3IRModel3wm_wrap.cxx -lstdc++ -lm -lluajit 
    