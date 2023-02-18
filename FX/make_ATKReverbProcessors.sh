    swig -lua -c++ -I../include -I../include ATKReverbProcessors.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKReverbProcessors.so ATKReverbProcessors_wrap.cxx -lstdc++ -lm -lluajit 
    