    swig -lua -c++ -Iinclude audiodsp_signal_mixer.swg
    gcc -Wfatal-errors -std=c++17 -Iinclude -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o audiodsp_signal_mixer.so audiodsp_signal_mixer_wrap.cxx -lstdc++ -lm -lluajit 
    