    swig -lua -c++ -I../include -I../include FxDSPRMSEstimator.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FxDSPRMSEstimator.so FxDSPRMSEstimator_wrap.cxx -lstdc++ -lm -lluajit 
    