    swig -lua -c++ -Iinclude VADPWOscillators.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o VADPWOscillators.so VADPWOscillators_wrap.cxx -lstdc++ -lm -lluajit    
    