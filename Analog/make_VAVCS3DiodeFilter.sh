    swig -lua -c++ -I../include -I../include/Analog VAVCS3DiodeFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAVCS3DiodeFilter.so VAVCS3DiodeFilter_wrap.cxx -lstdc++ -lm -lluajit 
    