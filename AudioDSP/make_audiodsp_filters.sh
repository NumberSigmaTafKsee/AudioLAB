    swig -lua -c++ -Iinclude audiodsp_filters.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_filters.so audiodsp_filters_wrap.cxx -lstdc++ -lm -lluajit 
    