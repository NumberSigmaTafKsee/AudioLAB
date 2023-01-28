    swig -lua -c++ -Iinclude FxDSPLinkwitzReillyFilter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FxDSPLinkwitzReillyFilter.so FxDSPLinkwitzReillyFilter_wrap.cxx -lstdc++ -lm -lluajit    
    