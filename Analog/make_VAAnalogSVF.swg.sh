    swig -lua -c++ -I../include -I../include/Analog VAAnalogSVF.swg.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o VAAnalogSVF.swg.so VAAnalogSVF.swg_wrap.cxx -lstdc++ -lm -lluajit 
    