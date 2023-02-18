    swig -lua -c++ -Iinclude audiodsp_voices.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_voices.so audiodsp_voices_wrap.cxx -lstdc++ -lm -lluajit 
    