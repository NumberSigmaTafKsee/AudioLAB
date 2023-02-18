    swig -lua -c++ -I../include -I../include/Analog VABlitOscillators.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VABlitOscillators.so VABlitOscillators_wrap.cxx -lstdc++ -lm -lluajit 
    