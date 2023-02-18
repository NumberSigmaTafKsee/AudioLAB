    swig -lua -c++ -Iinclude audiodsp_dynamics.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_dynamics.so audiodsp_dynamics_wrap.cxx -lstdc++ -lm -lluajit 
    