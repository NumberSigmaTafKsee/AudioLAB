    swig -lua -c++ -Iinclude LiquidNeuron.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o LiquidNeuron.so LiquidNeuron_wrap.cxx -lstdc++ -lm -lluajit    
    