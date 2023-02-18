    swig -lua -c++ -I../include -I../include/Analog VAStateVariableFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAStateVariableFilter.so VAStateVariableFilter_wrap.cxx -lstdc++ -lm -lluajit 
    