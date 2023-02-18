    swig -lua -c++ -Iinclude audiodsp_lv2plugin.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_lv2plugin.so audiodsp_lv2plugin_wrap.cxx -lstdc++ -lm -lluajit 
    