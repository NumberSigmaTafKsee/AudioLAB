    swig -lua -c++ -Iinclude BodeShifter.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o BodeShifter.so BodeShifter_wrap.cxx -lstdc++ -lm -lluajit    
    