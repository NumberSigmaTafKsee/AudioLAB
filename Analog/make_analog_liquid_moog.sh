    swig -lua -c++ -Iinclude analog_liquid_moog.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_liquid_moog.so analog_liquid_moog_wrap.cxx -lstdc++ -lm -lluajit 
    