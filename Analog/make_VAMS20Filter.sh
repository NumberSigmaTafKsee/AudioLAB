    swig -lua -c++ -I../include -I../include/Analog VAMS20Filter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMS20Filter.so VAMS20Filter_wrap.cxx -lstdc++ -lm -lluajit 
    