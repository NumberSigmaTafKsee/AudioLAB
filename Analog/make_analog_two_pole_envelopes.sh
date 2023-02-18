    swig -lua -c++ -Iinclude analog_two_pole_envelopes.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_two_pole_envelopes.so analog_two_pole_envelopes_wrap.cxx -lstdc++ -lm -lluajit 
    