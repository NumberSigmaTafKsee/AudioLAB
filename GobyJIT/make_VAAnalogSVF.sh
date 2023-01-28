    swig -lua -c++ -Iinclude VAAnalogSVF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VAAnalogSVF.so VAAnalogSVF_wrap.cxx -lstdc++ -lm -lluajit    
    