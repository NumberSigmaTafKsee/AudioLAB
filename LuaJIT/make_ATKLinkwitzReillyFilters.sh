    swig -lua -c++ -Iinclude ATKLinkwitzReillyFilters.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKLinkwitzReillyFilters.so ATKLinkwitzReillyFilters_wrap.cxx -lstdc++ -lm -lluajit    
    