    swig -lua -c++ -Iinclude Mu45FilterCalc.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Mu45FilterCalc.so Mu45FilterCalc_wrap.cxx -lstdc++ -lm -lluajit    
    