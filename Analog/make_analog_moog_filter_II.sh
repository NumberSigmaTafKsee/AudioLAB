    swig -lua -c++ -Iinclude analog_moog_filter_II.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_moog_filter_II.so analog_moog_filter_II_wrap.cxx -lstdc++ -lm -lluajit 
    