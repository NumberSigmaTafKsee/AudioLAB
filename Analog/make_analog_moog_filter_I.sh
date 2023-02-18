    swig -lua -c++ -Iinclude analog_moog_filter_I.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_moog_filter_I.so analog_moog_filter_I_wrap.cxx -lstdc++ -lm -lluajit 
    