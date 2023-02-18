    swig -lua -c++ -Iinclude analog_vco.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_vco.so analog_vco_wrap.cxx -lstdc++ -lm -lluajit 
    