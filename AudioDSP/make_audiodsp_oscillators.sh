    swig -lua -c++ -Iinclude audiodsp_oscillators.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_oscillators.so audiodsp_oscillators_wrap.cxx -lstdc++ -lm -lluajit 
    