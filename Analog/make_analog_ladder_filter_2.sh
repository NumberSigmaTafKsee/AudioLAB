    swig -lua -c++ -Iinclude analog_ladder_filter_2.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_ladder_filter_2.so analog_ladder_filter_2_wrap.cxx -lstdc++ -lm -lluajit 
    