    swig -lua -c++ -I../include -I../include/Analog VASVFSmoother.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VASVFSmoother.so VASVFSmoother_wrap.cxx -lstdc++ -lm -lluajit 
    