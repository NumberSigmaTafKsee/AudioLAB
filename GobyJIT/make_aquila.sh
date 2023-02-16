swig -lua -c++ -Iinclude/aquila aquila.i
gcc -Iinclude/aquila -O2 -fPIC -march=native -mavx2 -shared -o aquila.so aquila_wrap.cxx fft4g.c libAquila.a -lstdc++ -lm -lluajit
