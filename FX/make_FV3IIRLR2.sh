    swig -lua -c++ -I../include -I../include FV3IIRLR2.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3IIRLR2.so FV3IIRLR2_wrap.cxx -lstdc++ -lm -lluajit 
    