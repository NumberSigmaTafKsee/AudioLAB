    swig -lua -c++ -I../include -I../include FXProcessors.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FXProcessors.so FXProcessors_wrap.cxx -lstdc++ -lm -lluajit 
    