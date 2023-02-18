    swig -lua -c++ -I../include -I../include FV3SComp.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3SComp.so FV3SComp_wrap.cxx -lstdc++ -lm -lluajit 
    