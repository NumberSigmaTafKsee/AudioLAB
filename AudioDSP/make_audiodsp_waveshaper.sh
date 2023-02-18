    swig -lua -c++ -Iinclude audiodsp_waveshaper.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_waveshaper.so audiodsp_waveshaper_wrap.cxx -lstdc++ -lm -lluajit 
    