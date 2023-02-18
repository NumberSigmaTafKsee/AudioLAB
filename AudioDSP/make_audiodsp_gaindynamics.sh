    swig -lua -c++ -Iinclude audiodsp_gaindynamics.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_gaindynamics.so audiodsp_gaindynamics_wrap.cxx -lstdc++ -lm -lluajit 
    