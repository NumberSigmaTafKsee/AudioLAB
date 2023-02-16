    swig -lua -c++ -Iinclude AnalogSVF.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o AnalogSVF.so AnalogSVF_wrap.cxx -lstdc++ -lm -lluajit    
    