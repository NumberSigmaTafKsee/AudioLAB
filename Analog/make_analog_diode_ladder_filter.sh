    swig -lua -c++ -Iinclude analog_diode_ladder_filter.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_diode_ladder_filter.so analog_diode_ladder_filter_wrap.cxx -lstdc++ -lm -lluajit 
    