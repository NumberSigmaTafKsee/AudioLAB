    swig -lua -c++ -Iinclude audiodsp_simple_compressor.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_simple_compressor.so audiodsp_simple_compressor_wrap.cxx -lstdc++ -lm -lluajit 
    