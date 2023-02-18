    swig -lua -c++ -Iinclude audiodsp_wavechild.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_wavechild.so audiodsp_wavechild_wrap.cxx -lstdc++ -lm -lluajit 
    