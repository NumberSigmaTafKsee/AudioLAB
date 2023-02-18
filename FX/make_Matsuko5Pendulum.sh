    swig -lua -c++ -I../include -I../include Matsuko5Pendulum.swg
    gcc -Wfatal-errors -std=c++17 -I.. -I../include -O2 -fPIC -march=native -mavx2 -fopenmp -pthread -shared -o Matsuko5Pendulum.so Matsuko5Pendulum_wrap.cxx -lstdc++ -lm -lluajit 
    