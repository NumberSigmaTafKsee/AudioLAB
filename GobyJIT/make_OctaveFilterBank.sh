    swig -lua -c++ -Iinclude OctaveFilterBank.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o OctaveFilterBank.so OctaveFilterBank_wrap.cxx -lstdc++ -lm -lluajit    
    