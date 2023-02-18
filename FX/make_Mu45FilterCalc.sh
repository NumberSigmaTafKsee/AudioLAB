    swig -lua -c++ -I../include -I../include Mu45FilterCalc.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Mu45FilterCalc.so Mu45FilterCalc_wrap.cxx -lstdc++ -lm -lluajit 
    