    swig -lua -c++ -I../include -I../include PeakHoldCascade2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o PeakHoldCascade2.so PeakHoldCascade2_wrap.cxx -lstdc++ -lm -lluajit 
    