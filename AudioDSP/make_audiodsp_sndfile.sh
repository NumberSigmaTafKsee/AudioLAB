    swig -lua -c++ -Iinclude audiodsp_sndfile.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_sndfile.so audiodsp_sndfile_wrap.cxx -lstdc++ -lm -lluajit 
    