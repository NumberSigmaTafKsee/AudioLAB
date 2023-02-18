    swig -lua -c++ -Iinclude analog_diode_clipper.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_diode_clipper.so analog_diode_clipper_wrap.cxx -lstdc++ -lm -lluajit 
    