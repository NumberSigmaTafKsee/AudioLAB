    swig -lua -c++ -I../include -I../include/Analog VAVCS3Filter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAVCS3Filter.so VAVCS3Filter_wrap.cxx -lstdc++ -lm -lluajit 
    