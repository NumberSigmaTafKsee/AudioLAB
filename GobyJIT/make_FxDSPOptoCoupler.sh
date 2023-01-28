    swig -lua -c++ -Iinclude FxDSPOptoCoupler.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPOptoCoupler.so FxDSPOptoCoupler_wrap.cxx -lstdc++ -lm -lluajit    
    