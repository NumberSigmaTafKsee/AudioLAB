    swig -lua -c++ -Iinclude audiodsp_amplifiers.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_amplifiers.so audiodsp_amplifiers_wrap.cxx -lstdc++ -lm -lluajit 
    