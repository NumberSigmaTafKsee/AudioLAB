swig -lua -c++ -Ispuce spuce.i
gcc -fmax-errors=1 -I. -Ispuce -O2 -fPIC -march=native -mavx2 -shared -o spuce.so spuce_wrap.cxx libspuce.a -lstdc++ -lm -lluajit
