    swig -lua -c++ -Iinclude audiodsp_doom_compressor.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_doom_compressor.so audiodsp_doom_compressor_wrap.cxx -lstdc++ -lm -lluajit 
    