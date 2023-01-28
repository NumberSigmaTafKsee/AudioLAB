    swig -lua -c++ -Iinclude FXProcessor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o FXProcessor.so FXProcessor_wrap.cxx -lstdc++ -lm -lluajit    
    