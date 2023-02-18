    swig -lua -c++ -Iinclude analog_blit_oscillators.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_blit_oscillators.so analog_blit_oscillators_wrap.cxx -lstdc++ -lm -lluajit 
    