    swig -lua -c++ -Iinclude audiodsp_wave_table.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_wave_table.so audiodsp_wave_table_wrap.cxx -lstdc++ -lm -lluajit 
    