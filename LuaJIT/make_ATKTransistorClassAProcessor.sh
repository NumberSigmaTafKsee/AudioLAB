    swig -lua -c++ -Iinclude ATKTransistorClassAProcessor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKTransistorClassAProcessor.so ATKTransistorClassAProcessor_wrap.cxx -lstdc++ -lm -lluajit    
    