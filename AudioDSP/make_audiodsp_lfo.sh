    swig -lua -c++ -Iinclude audiodsp_lfo.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_lfo.so audiodsp_lfo_wrap.cxx -lstdc++ -lm -lluajit 
    