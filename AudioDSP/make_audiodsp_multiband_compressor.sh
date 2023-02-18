    swig -lua -c++ -Iinclude audiodsp_multiband_compressor.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_multiband_compressor.so audiodsp_multiband_compressor_wrap.cxx -lstdc++ -lm -lluajit 
    