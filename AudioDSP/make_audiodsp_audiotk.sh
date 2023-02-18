    swig -lua -c++ -Iinclude audiodsp_audiotk.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_audiotk.so audiodsp_audiotk_wrap.cxx -lstdc++ -lm -lluajit 
    