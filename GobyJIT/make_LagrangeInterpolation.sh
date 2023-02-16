    swig -lua -c++ -Iinclude LagrangeInterpolation.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o LagrangeInterpolation.so LagrangeInterpolation_wrap.cxx -lstdc++ -lm -lluajit    
    