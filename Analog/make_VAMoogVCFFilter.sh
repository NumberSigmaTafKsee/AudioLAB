    swig -lua -c++ -I../include -I../include/Analog VAMoogVCFFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogVCFFilter.so VAMoogVCFFilter_wrap.cxx -lstdc++ -lm -lluajit 
    