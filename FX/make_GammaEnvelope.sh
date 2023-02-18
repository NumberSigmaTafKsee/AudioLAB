    swig -lua -c++ -I../include -I../include GammaEnvelope.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o GammaEnvelope.so GammaEnvelope_wrap.cxx -lstdc++ -lm -lluajit 
    