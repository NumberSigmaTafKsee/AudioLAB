    swig -lua -c++ -I../include -I../include RackFXEqualizer.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o RackFXEqualizer.so RackFXEqualizer_wrap.cxx -lstdc++ -lm -lluajit 
    