    swig -lua -c++ -I../include -I../include TinyTricksOscillators.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o TinyTricksOscillators.so TinyTricksOscillators_wrap.cxx -lstdc++ -lm -lluajit 
    