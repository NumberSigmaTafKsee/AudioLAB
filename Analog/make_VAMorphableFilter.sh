    swig -lua -c++ -I../include -I../include/Analog VAMorphableFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMorphableFilter.so VAMorphableFilter_wrap.cxx -lstdc++ -lm -lluajit 
    