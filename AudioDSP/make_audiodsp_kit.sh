    swig -lua -c++ -Iinclude audiodsp_kit.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_kit.so audiodsp_kit_wrap.cxx -lstdc++ -lm -lluajit 
    