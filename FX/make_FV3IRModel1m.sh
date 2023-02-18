    swig -lua -c++ -I../include -I../include FV3IRModel1m.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IRModel1m.so FV3IRModel1m_wrap.cxx -lstdc++ -lm -lluajit 
    