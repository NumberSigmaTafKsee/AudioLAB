    swig -lua -c++ -I../include -I../include FXCE2Chorus.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FXCE2Chorus.so FXCE2Chorus_wrap.cxx -lstdc++ -lm -lluajit 
    