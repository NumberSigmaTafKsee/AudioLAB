    swig -lua -c++ -Iinclude analog_sst_waveshapers.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_sst_waveshapers.so analog_sst_waveshapers_wrap.cxx -lstdc++ -lm -lluajit 
    