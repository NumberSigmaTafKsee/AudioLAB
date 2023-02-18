    swig -lua -c++ -I../include -I../include ATKTriode2Processor.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKTriode2Processor.so ATKTriode2Processor_wrap.cxx -lstdc++ -lm -lluajit 
    