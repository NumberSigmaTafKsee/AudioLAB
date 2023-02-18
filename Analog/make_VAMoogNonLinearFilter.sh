    swig -lua -c++ -I../include -I../include/Analog VAMoogNonLinearFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogNonLinearFilter.so VAMoogNonLinearFilter_wrap.cxx -lstdc++ -lm -lluajit 
    