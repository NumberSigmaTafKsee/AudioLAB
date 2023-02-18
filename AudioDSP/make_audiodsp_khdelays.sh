    swig -lua -c++ -Iinclude audiodsp_khdelays.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_khdelays.so audiodsp_khdelays_wrap.cxx -lstdc++ -lm -lluajit 
    