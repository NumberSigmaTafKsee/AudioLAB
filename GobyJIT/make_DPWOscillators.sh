    swig -lua -c++ -Iinclude DPWOscillators.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o DPWOscillators.so DPWOscillators_wrap.cxx -lstdc++ -lm -lluajit    
    