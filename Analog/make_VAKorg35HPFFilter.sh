    swig -lua -c++ -I../include -I../include/Analog VAKorg35HPFFilter.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAKorg35HPFFilter.so VAKorg35HPFFilter_wrap.cxx -lstdc++ -lm -lluajit 
    