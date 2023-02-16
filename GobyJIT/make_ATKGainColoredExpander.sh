    swig -lua -c++ -Iinclude ATKGainColoredExpander.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKGainColoredExpander.so ATKGainColoredExpander_wrap.cxx -lstdc++ -lm -lluajit    
    