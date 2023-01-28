    swig -lua -c++ -Iinclude FxDSPRMSEstimator.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPRMSEstimator.so FxDSPRMSEstimator_wrap.cxx -lstdc++ -lm -lluajit    
    