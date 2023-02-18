    swig -lua -c++ -Iinclude analog_kit.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_kit.so analog_kit_wrap.cxx -lstdc++ -lm -lluajit 
    