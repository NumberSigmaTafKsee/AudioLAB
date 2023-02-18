    swig -lua -c++ -Iinclude audiodsp_compressor_base.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_compressor_base.so audiodsp_compressor_base_wrap.cxx -lstdc++ -lm -lluajit 
    