    swig -lua -c++ -Iinclude AudioDSP_ModDelay.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AudioDSP_ModDelay.so AudioDSP_ModDelay_wrap.cxx -lstdc++ -lm -lluajit    
    