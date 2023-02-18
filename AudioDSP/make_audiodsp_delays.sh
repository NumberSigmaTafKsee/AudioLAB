    swig -lua -c++ -Iinclude audiodsp_delays.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_delays.so audiodsp_delays_wrap.cxx -lstdc++ -lm -lluajit 
    