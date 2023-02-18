    swig -lua -c++ -Iinclude audiodsp_ctagdrc.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_ctagdrc.so audiodsp_ctagdrc_wrap.cxx -lstdc++ -lm -lluajit 
    