    swig -lua -c++ -I../include -I../include/Analog VAMoogLadderOberheim.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogLadderOberheim.so VAMoogLadderOberheim_wrap.cxx -lstdc++ -lm -lluajit 
    