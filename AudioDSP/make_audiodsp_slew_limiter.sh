    swig -lua -c++ -Iinclude audiodsp_slew_limiter.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_slew_limiter.so audiodsp_slew_limiter_wrap.cxx -lstdc++ -lm -lluajit 
    