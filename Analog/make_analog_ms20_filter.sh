    swig -lua -c++ -Iinclude analog_ms20_filter.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_ms20_filter.so analog_ms20_filter_wrap.cxx -lstdc++ -lm -lluajit 
    