    swig -lua -c++ -Iinclude analog_wdf_passive_lpf.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_wdf_passive_lpf.so analog_wdf_passive_lpf_wrap.cxx -lstdc++ -lm -lluajit 
    