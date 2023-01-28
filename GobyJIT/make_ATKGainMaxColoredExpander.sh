    swig -lua -c++ -Iinclude ATKGainMaxColoredExpander.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKGainMaxColoredExpander.so ATKGainMaxColoredExpander_wrap.cxx -lstdc++ -lm -lluajit    
    