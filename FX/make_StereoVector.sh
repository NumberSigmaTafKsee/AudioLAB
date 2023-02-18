    swig -lua -c++ -I../include -I../include StereoVector.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o StereoVector.so StereoVector_wrap.cxx -lstdc++ -lm -lluajit 
    