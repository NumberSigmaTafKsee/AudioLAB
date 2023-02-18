    swig -lua -c++ -Iinclude audiodsp_ccaubio.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_ccaubio.so audiodsp_ccaubio_wrap.cxx -lstdc++ -lm -lluajit 
    