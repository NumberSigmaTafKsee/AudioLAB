    swig -lua -c++ -I../include -I../include FV3RevModel.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o FV3RevModel.so FV3RevModel_wrap.cxx -lstdc++ -lm -lluajit 
    