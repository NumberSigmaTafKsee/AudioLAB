    swig -lua -c++ -I../include -I../include NOBSNonlinearOscillator.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o NOBSNonlinearOscillator.so NOBSNonlinearOscillator_wrap.cxx -lstdc++ -lm -lluajit 
    