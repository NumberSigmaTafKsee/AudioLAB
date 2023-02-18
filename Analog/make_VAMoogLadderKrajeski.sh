    swig -lua -c++ -I../include -I../include/Analog VAMoogLadderKrajeski.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogLadderKrajeski.so VAMoogLadderKrajeski_wrap.cxx -lstdc++ -lm -lluajit 
    