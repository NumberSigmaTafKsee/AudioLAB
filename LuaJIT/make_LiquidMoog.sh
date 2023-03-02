    swig -lua -c++ -Iinclude LiquidMoog.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o LiquidMoog.so LiquidMoog_wrap.cxx -lstdc++ -lm -lluajit    
    