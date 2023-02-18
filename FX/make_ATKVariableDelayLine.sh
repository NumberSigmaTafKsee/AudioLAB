    swig -lua -c++ -I../include -I../include ATKVariableDelayLine.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o ATKVariableDelayLine.so ATKVariableDelayLine_wrap.cxx -lstdc++ -lm -lluajit 
    