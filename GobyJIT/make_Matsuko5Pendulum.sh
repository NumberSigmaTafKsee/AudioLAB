    swig -lua -c++ -Iinclude Matsuko5Pendulum.i
    gcc -std=c++17 -I. -I.. -Iinclude -O2 -fPIC -march=native -mavx2 -shared -o Matsuko5Pendulum.so Matsuko5Pendulum_wrap.cxx -lstdc++ -lm -lluajit    
    