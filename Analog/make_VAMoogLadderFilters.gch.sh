    swig -lua -c++ -I../include -I../include/Analog VAMoogLadderFilters.gch.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogLadderFilters.gch.so VAMoogLadderFilters.gch_wrap.cxx -lstdc++ -lm -lluajit 
    