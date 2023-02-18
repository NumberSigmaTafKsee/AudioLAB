    swig -lua -c++ -Iinclude analog_dpw_saw_oscillator.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_dpw_saw_oscillator.so analog_dpw_saw_oscillator_wrap.cxx -lstdc++ -lm -lluajit 
    