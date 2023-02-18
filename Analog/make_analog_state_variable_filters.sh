    swig -lua -c++ -Iinclude analog_state_variable_filters.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_state_variable_filters.so analog_state_variable_filters_wrap.cxx -lstdc++ -lm -lluajit 
    