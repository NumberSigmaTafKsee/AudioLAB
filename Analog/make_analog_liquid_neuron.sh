    swig -lua -c++ -Iinclude analog_liquid_neuron.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_liquid_neuron.so analog_liquid_neuron_wrap.cxx -lstdc++ -lm -lluajit 
    