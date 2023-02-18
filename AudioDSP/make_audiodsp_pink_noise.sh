    swig -lua -c++ -Iinclude audiodsp_pink_noise.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_pink_noise.so audiodsp_pink_noise_wrap.cxx -lstdc++ -lm -lluajit 
    