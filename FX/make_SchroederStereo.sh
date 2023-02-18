    swig -lua -c++ -I../include -I../include SchroederStereo.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o SchroederStereo.so SchroederStereo_wrap.cxx -lstdc++ -lm -lluajit 
    