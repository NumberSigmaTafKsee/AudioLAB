    swig -lua -c++ -I../include -I../include/Analog VAMoogLadders.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogLadders.so VAMoogLadders_wrap.cxx -lstdc++ -lm -lluajit 
    