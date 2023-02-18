    swig -lua -c++ -Iinclude audiodsp_compressor.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_compressor.so audiodsp_compressor_wrap.cxx -lstdc++ -lm -lluajit 
    