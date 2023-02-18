    swig -lua -c++ -Iinclude analog_svf_smooth_filter.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_svf_smooth_filter.so analog_svf_smooth_filter_wrap.cxx -lstdc++ -lm -lluajit 
    