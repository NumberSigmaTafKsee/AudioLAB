    swig -lua -c++ -I../include -I../include ATKGainColoredExpander.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKGainColoredExpander.so ATKGainColoredExpander_wrap.cxx -lstdc++ -lm -lluajit 
    