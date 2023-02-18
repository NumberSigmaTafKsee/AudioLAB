    swig -lua -c++ -I../include -I../include LiquidMoog.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o LiquidMoog.so LiquidMoog_wrap.cxx -lstdc++ -lm -lluajit 
    