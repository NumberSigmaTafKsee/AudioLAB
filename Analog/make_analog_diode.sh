    swig -lua -c++ -Iinclude analog_diode.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_diode.so analog_diode_wrap.cxx -lstdc++ -lm -lluajit 
    