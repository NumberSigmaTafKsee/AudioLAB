gcc -O2 -fPIC -march=native -mavx2 -c *.cpp
ar -rcv -o libfftconvolver.a *.o
