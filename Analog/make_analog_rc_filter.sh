    swig -lua -c++ -Iinclude analog_rc_filter.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_rc_filter.so analog_rc_filter_wrap.cxx -lstdc++ -lm -lluajit 
    