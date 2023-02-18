    swig -lua -c++ -Iinclude audiodsp_limiter.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_limiter.so audiodsp_limiter_wrap.cxx -lstdc++ -lm -lluajit 
    