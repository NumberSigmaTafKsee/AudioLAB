    swig -lua -c++ -Iinclude analog_slew_limiter.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_slew_limiter.so analog_slew_limiter_wrap.cxx -lstdc++ -lm -lluajit 
    