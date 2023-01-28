swig -lua -c++ -Iinclude Std.i
gcc -O2 -Iinclude/Std -fPIC -march=native -mavx2 -c *.cpp
gcc -pthread -std=c++17 -fmax-errors=1 -Iinclude -Iinclude/Std -O2 -fPIC -shared -o Std.so Std_wrap.cxx include/gnuplot_i.c *.o -lstdc++ -lstdc++fs -lm -lluajit -L/usr/local/lib -ljsoncpp -lrt 
