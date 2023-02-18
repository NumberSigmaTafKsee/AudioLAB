    swig -lua -c++ -I../include -I../include/Analog VASstWaveshaper.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VASstWaveshaper.so VASstWaveshaper_wrap.cxx -lstdc++ -lm -lluajit 
    