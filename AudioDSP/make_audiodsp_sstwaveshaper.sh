    swig -lua -c++ -Iinclude audiodsp_sstwaveshaper.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_sstwaveshaper.so audiodsp_sstwaveshaper_wrap.cxx -lstdc++ -lm -lluajit 
    