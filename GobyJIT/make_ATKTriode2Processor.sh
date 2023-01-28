    swig -lua -c++ -Iinclude ATKTriode2Processor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKTriode2Processor.so ATKTriode2Processor_wrap.cxx -lstdc++ -lm -lluajit    
    