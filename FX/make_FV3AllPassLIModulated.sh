    swig -lua -c++ -I../include -I../include FV3AllPassLIModulated.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3AllPassLIModulated.so FV3AllPassLIModulated_wrap.cxx -lstdc++ -lm -lluajit 
    