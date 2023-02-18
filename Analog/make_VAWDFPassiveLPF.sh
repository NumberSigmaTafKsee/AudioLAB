    swig -lua -c++ -I../include -I../include/Analog VAWDFPassiveLPF.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAWDFPassiveLPF.so VAWDFPassiveLPF_wrap.cxx -lstdc++ -lm -lluajit 
    