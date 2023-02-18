    swig -lua -c++ -I../include -I../include FXFreeVerb3.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FXFreeVerb3.so FXFreeVerb3_wrap.cxx -lstdc++ -lm -lluajit 
    