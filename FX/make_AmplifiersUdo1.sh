    swig -lua -c++ -I../include -I../include AmplifiersUdo1.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o AmplifiersUdo1.so AmplifiersUdo1_wrap.cxx -lstdc++ -lm -lluajit 
    