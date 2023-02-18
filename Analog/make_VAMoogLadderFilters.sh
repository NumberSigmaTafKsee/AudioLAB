    swig -lua -c++ -I../include -I../include/Analog VAMoogLadderFilters.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogLadderFilters.so VAMoogLadderFilters_wrap.cxx -lstdc++ -lm -lluajit 
    