    swig -lua -c++ -Iinclude audiodsp_vafilters.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_vafilters.so audiodsp_vafilters_wrap.cxx -lstdc++ -lm -lluajit 
    