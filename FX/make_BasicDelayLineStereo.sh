    swig -lua -c++ -I../include -I../include BasicDelayLineStereo.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o BasicDelayLineStereo.so BasicDelayLineStereo_wrap.cxx -lstdc++ -lm -lluajit 
    