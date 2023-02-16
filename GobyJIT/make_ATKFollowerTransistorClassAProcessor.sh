    swig -lua -c++ -Iinclude ATKFollowerTransistorClassAProcessor.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o ATKFollowerTransistorClassAProcessor.so ATKFollowerTransistorClassAProcessor_wrap.cxx -lstdc++ -lm -lluajit    
    