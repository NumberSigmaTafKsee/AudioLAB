    swig -lua -c++ -I../include -I../include FV3IRModel2z.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IRModel2z.so FV3IRModel2z_wrap.cxx -lstdc++ -lm -lluajit 
    