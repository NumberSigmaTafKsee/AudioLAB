    swig -lua -c++ -Iinclude analog_svf_smoother.swg
    gcc -Wfatal-errors -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o analog_svf_smoother.so analog_svf_smoother_wrap.cxx -lstdc++ -lm -lluajit 
    