    swig -lua -c++ -Iinclude analog_runge_kutta_filter.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_runge_kutta_filter.so analog_runge_kutta_filter_wrap.cxx -lstdc++ -lm -lluajit 
    