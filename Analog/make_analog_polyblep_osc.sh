    swig -lua -c++ -Iinclude analog_polyblep_osc.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_polyblep_osc.so analog_polyblep_osc_wrap.cxx -lstdc++ -lm -lluajit 
    