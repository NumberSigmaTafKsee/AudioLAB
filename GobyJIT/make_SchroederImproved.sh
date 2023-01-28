    swig -lua -c++ -Iinclude SchroederImproved.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o SchroederImproved.so SchroederImproved_wrap.cxx -lstdc++ -lm -lluajit    
    