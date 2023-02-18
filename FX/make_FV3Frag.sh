    swig -lua -c++ -I../include -I../include FV3Frag.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3Frag.so FV3Frag_wrap.cxx -lstdc++ -lm -lluajit 
    