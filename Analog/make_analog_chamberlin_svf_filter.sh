    swig -lua -c++ -Iinclude analog_chamberlin_svf_filter.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_chamberlin_svf_filter.so analog_chamberlin_svf_filter_wrap.cxx -lstdc++ -lm -lluajit 
    