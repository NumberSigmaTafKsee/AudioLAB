    swig -lua -c++ -I../include -I../include/Analog VAMoogNonLinearFilter2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAMoogNonLinearFilter2.so VAMoogNonLinearFilter2_wrap.cxx -lstdc++ -lm -lluajit 
    