    swig -lua -c++ -Iinclude audiodsp_sstfilters.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_sstfilters.so audiodsp_sstfilters_wrap.cxx -lstdc++ -lm -lluajit 
    