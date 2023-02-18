    swig -lua -c++ -Iinclude audiodsp_low_frequency_oscillator.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_low_frequency_oscillator.so audiodsp_low_frequency_oscillator_wrap.cxx -lstdc++ -lm -lluajit 
    