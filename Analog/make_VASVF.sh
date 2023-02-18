    swig -lua -c++ -I../include -I../include/Analog VASVF.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VASVF.so VASVF_wrap.cxx -lstdc++ -lm -lluajit 
    