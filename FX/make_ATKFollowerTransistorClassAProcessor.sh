    swig -lua -c++ -I../include -I../include ATKFollowerTransistorClassAProcessor.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKFollowerTransistorClassAProcessor.so ATKFollowerTransistorClassAProcessor_wrap.cxx -lstdc++ -lm -lluajit 
    