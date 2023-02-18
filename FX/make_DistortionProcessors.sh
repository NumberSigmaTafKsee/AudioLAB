    swig -lua -c++ -I../include -I../include DistortionProcessors.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o DistortionProcessors.so DistortionProcessors_wrap.cxx -lstdc++ -lm -lluajit 
    